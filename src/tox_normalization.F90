#include "macros.h"

!> Module with normalization routines for tensor omics.
module tox_normalization
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, set_err, ERR_EMPTY_INPUT, ERR_DIVISION_BY_ZERO, ERR_INVALID_INPUT, is_err, validate_dimension_size, validate_in_range_real, ERR_ALLOC_FAIL, validate_all_in_range_int
    use f42_utils, only: norm, is_close, logx, mean, std_dev
    use tox_loess, only: loess_alloc

#define CM_LOESS_SPAN_DEFAULT 0.7_real64
#define CM_LOESS_DEGREE_DEFAULT 2_int32

contains

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Normalizes an input vector to unit length in-place
    pure subroutine normalize_unit_length(vector, n_dims, ierr)
        integer(int32), intent(in) :: n_dims
            !! number of elements in `vector`
        real(real64), dimension(n_dims), intent(inout) :: vector
            !! Vector that will be normalized to unit length
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_dim
        real(real64) :: vector_norm

        call set_ok(ierr)

        call validate_dimension_size(n_dims, ierr)
        if (is_err(ierr)) return

        vector_norm = norm(vector)
        if (is_close(vector_norm, 0.0_real64)) then
            call set_err(ierr, ERR_DIVISION_BY_ZERO)
            return
        end if

        ! check for nan, inf
        call validate_in_range_real(vector_norm, ierr)
        if (is_err(ierr)) return

        do concurrent (i_dim = 1:n_dims) shared(vector, vector_norm)
            vector(i_dim) = vector(i_dim)/vector_norm
        end do
    end subroutine normalize_unit_length

    !> AUTHOR_VIVIAN_BASS
    !| Complete normalization pipeline for gene expression data.
    !| Performs: std dev normalization, quantile normalization, replicate averaging, log2(x+1) transformation.
    !| Final result is in log_transformed_expr. If fold change is needed, call calc_fchange separately.
    subroutine normalization_pipeline_alloc(n_genes, n_replicates, expr, log_transformed_expr, reps_per_tissue, n_tissues, span, degree, use_quantile, ierr)

        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
        integer(int32), dimension(n_tissues), intent(in) :: reps_per_tissue
            !! Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
            !! e.g. `[2,3]` means `5` total replicates per gene, with the `expr(1:2, i_gene)` related to the first tissue and `expr(3:, i_gene)` related to the second one.
        real(real64), dimension(n_tissues, n_genes), intent(out), target :: log_transformed_expr
            !! Log-transformed grouped `expr`

        real(real64), intent(in), optional :: span
            !! LOESS span parameter, default: `CM_LOESS_SPAN_DEFAULT`
        integer(int32), intent(in), optional :: degree
            !! LOESS degree parameter, default: `CM_LOESS_DEGREE_DEFAULT`
        logical, intent(in), optional :: use_quantile
            !! Use quantile normalization, default: `.false.`
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Local variables
        logical :: actual_use_quantile

        real(real64), dimension(:), allocatable, target :: tmp_loess_y, tmp_yhat_global
        real(real64), dimension(:, :), allocatable :: expr_copy
        integer(int32), dimension(:), allocatable, target :: tmp_indices_used

        real(real64), dimension(:, :), pointer :: log_transformed_expr_transposed_view
        real(real64), dimension(:), pointer :: tmp_loess_x_ptr, tmp_loess_y_ptr, tmp_yhat_global_ptr, tmp_genes_row_ptr, tmp_rank_means_ptr
        integer(int32), dimension(:), pointer :: tmp_indices_used_ptr, tmp_perm_ptr

        ! Error handling
        call set_ok(ierr)

        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_replicates, ierr)
        call validate_dimension_size(n_tissues, ierr)

        if (is_err(ierr)) return

        M_DEFAULT_VAL(use_quantile, actual_use_quantile, .false.)

        log_transformed_expr_transposed_view(1:n_genes, 1:n_tissues) => log_transformed_expr
        tmp_loess_x_ptr => log_transformed_expr_transposed_view(:, 1)
        select case (n_tissues)
            case (1)
                M_ALLOCATE(tmp_loess_y(n_genes))
                tmp_loess_y_ptr => tmp_loess_y
                M_ALLOCATE(tmp_yhat_global(n_genes))
                tmp_yhat_global_ptr => tmp_yhat_global
            case (2)
                tmp_loess_y_ptr => log_transformed_expr_transposed_view(:, 2)
                M_ALLOCATE(tmp_yhat_global(n_genes))
                tmp_yhat_global_ptr => tmp_yhat_global
            case (3:)
                tmp_loess_y_ptr => log_transformed_expr_transposed_view(:, 2)
                tmp_yhat_global_ptr => log_transformed_expr_transposed_view(:, 3)
        end select

        M_ALLOCATE(tmp_indices_used(n_genes))
        tmp_indices_used_ptr => tmp_indices_used

        M_ALLOCATE(expr_copy(n_replicates, n_genes))
        expr_copy = expr

        ! Step 1: Normalize per-gene by std dev
        call normalize_by_std_dev_inplace_helper(n_genes, n_replicates, expr_copy, tmp_loess_x_ptr, tmp_loess_y_ptr, tmp_indices_used_ptr, tmp_yhat_global_ptr, span, degree, ierr)
        if (is_err(ierr)) return

        ! Step 2: Quantile normalization (conditional)
        if (actual_use_quantile) then
            tmp_perm_ptr => tmp_indices_used_ptr
            tmp_genes_row_ptr => tmp_loess_x_ptr
            tmp_rank_means_ptr => tmp_loess_y_ptr
            call quantile_normalization_inplace_helper(n_genes, n_replicates, expr_copy, tmp_genes_row_ptr, tmp_rank_means_ptr, tmp_perm_ptr)
        end if

        call calc_tiss_avg_helper(n_genes, n_tissues, reps_per_tissue, expr_copy, log_transformed_expr)
        if (is_err(ierr)) return

        ! Step 4: Log2(x+1) transformation
        call log2_transformation_inplace_helper(n_genes, n_tissues, log_transformed_expr, ierr)
        if (is_err(ierr)) return
    end subroutine normalization_pipeline_alloc

    !> AUTHOR_VIVIAN_BASS
    !| Normalizes each gene's expression vector using LOESS-stabilized standard deviation.
    !| This procedure applies a global stabilization based on the relationship between
    !| gene-wise mean expression and empirical standard deviation.
    subroutine normalize_by_std_dev_alloc(n_genes, n_replicates, expr, normalized_expr, span, degree, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
        real(real64), dimension(n_replicates, n_genes), intent(out) :: normalized_expr
            !! Normalized `expr`
        real(real64), intent(in), optional :: span
            !! LOESS span parameter, default: `CM_LOESS_SPAN_DEFAULT`
        integer(int32), intent(in), optional :: degree
            !! LOESS degree parameter, default: `CM_LOESS_DEGREE_DEFAULT`
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Buffers for LOESS fitting (preallocated to avoid internal allocations)
        real(real64), dimension(:), allocatable :: tmp_loess_x, tmp_loess_y, tmp_yhat_global
        integer(int32), dimension(:), allocatable :: tmp_indices_used

        M_ALLOCATE(tmp_loess_x(n_genes))
        M_ALLOCATE(tmp_loess_y(n_genes))
        M_ALLOCATE(tmp_yhat_global(n_genes))
        M_ALLOCATE(tmp_indices_used(n_genes))

        normalized_expr = expr
        call normalize_by_std_dev_inplace_helper(n_genes, n_replicates, normalized_expr, &
                                    tmp_loess_x, tmp_loess_y, tmp_indices_used, tmp_yhat_global, &
                                    span, degree, ierr)
    end subroutine normalize_by_std_dev_alloc

    !> AUTHOR_VIVIAN_BASS
    !| Normalizes each gene's expression vector using LOESS-stabilized standard deviation.
    !| This procedure applies a global stabilization based on the relationship between
    !| gene-wise mean expression and empirical standard deviation.
    subroutine normalize_by_std_dev_inplace_helper(n_genes, n_replicates, expr, &
                                    tmp_loess_x, tmp_loess_y, tmp_indices_used, tmp_yhat_global, &
                                    span, degree, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(inout) :: expr
            !! Gene Expression matrix

        ! Buffers for LOESS fitting (preallocated to avoid internal allocations)
        real(real64), dimension(n_genes), intent(out) :: tmp_loess_x
            !! Mean values (X-axis for LOESS)
        real(real64), dimension(n_genes), intent(out) :: tmp_loess_y
            !! Empirical standard deviation values (Y-axis for LOESS)
        integer(int32), dimension(n_genes), intent(out) :: tmp_indices_used
            !! Mapping back to gene indices
        real(real64), dimension(n_genes), intent(out) :: tmp_yhat_global
            !! Fitted standard deviation values (LOESS predictions)

        real(real64), intent(in), optional :: span
            !! LOESS span parameter, default: `CM_LOESS_SPAN_DEFAULT`
        integer(int32), intent(in), optional :: degree
            !! LOESS degree parameter, default: `CM_LOESS_DEGREE_DEFAULT`
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Local variables
        integer(int32) :: i_gene, i_valid, i_tissue, n_valid, gene_idx, actual_degree
        real(real64) :: mean_val, fitted_sd, actual_span

        M_DEFAULT_VAL(span, actual_span, CM_LOESS_SPAN_DEFAULT)
        M_DEFAULT_VAL(degree, actual_degree, CM_LOESS_DEGREE_DEFAULT)

        ! Initialize error code and output arrays
        call set_ok(ierr)
        n_valid = 0
        tmp_yhat_global = 0.0_real64

        ! Step 1: stats per gene
        do i_gene = 1, n_genes
            mean_val = mean(expr(:, i_gene))
            tmp_loess_x(i_gene) = mean_val
            tmp_loess_y(i_gene) = std_dev(expr(:, i_gene))

            if (is_close(tmp_loess_y(i_gene), 0.0_real64)) cycle

            n_valid = n_valid + 1
            tmp_loess_x(n_valid) = tmp_loess_x(i_gene)
            tmp_loess_y(n_valid) = tmp_loess_y(i_gene)
            tmp_indices_used(n_valid) = i_gene
        end do

        if (n_valid < 5) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        call loess_alloc(x=tmp_loess_x(1:n_valid), y=tmp_loess_y(1:n_valid), &
                         span=actual_span, degree=actual_degree, fitted_values=tmp_yhat_global(1:n_valid), &
                         mode=1, n_iters=3, ierr=ierr)
        if (is_err(ierr)) return

        ! Step 3: apply normalization
        do concurrent (i_valid = 1:n_valid) local(fitted_sd, gene_idx) shared(tmp_yhat_global, tmp_loess_y, tmp_indices_used, expr)
            fitted_sd = tmp_yhat_global(i_valid)
            if (is_close(fitted_sd, 0.0_real64)) fitted_sd = tmp_loess_y(i_valid)

            gene_idx = tmp_indices_used(i_valid)
            do concurrent (i_tissue = 1:n_replicates) shared (expr, gene_idx, fitted_sd)
                expr(i_tissue, gene_idx) = expr(i_tissue, gene_idx) / fitted_sd
            end do
        end do
    end subroutine normalize_by_std_dev_inplace_helper

    !> AUTHOR_VIVIAN_BASS
    !| Normalizes each gene's expression vector using `sqrt(mean(x^2))`
    !| across tissues (not classical standard deviation).
    pure subroutine root_mean_sq_normalization(n_genes, n_replicates, expr, normalized_expr, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
        real(real64), dimension(n_replicates, n_genes), intent(out) :: normalized_expr
            !! Normalized `expr`
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Error handling
        call set_ok(ierr)

        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_replicates, ierr)

        if (is_err(ierr)) return

        normalized_expr = expr
        call root_mean_sq_normalization_inplace_helper(n_genes, n_replicates, normalized_expr)
    end subroutine root_mean_sq_normalization

    !> AUTHOR_VIVIAN_BASS
    !| Normalizes each gene's expression vector using `sqrt(mean(x^2))`
    !| across tissues (not classical standard deviation).
    pure subroutine root_mean_sq_normalization_inplace_helper(n_genes, n_replicates, expr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(inout) :: expr
            !! Gene Expression matrix

        ! Local variables
        integer(int32) :: i_gene, i_tissue
        real(real64) :: std_dev, temp_sum

        ! Loop over each gene
        do concurrent (i_gene = 1:n_genes) local(temp_sum, std_dev) shared(n_replicates, expr)
            temp_sum = 0.0_real64
            do concurrent (i_tissue = 1:n_replicates) shared(expr, i_gene) reduce(+:temp_sum)
                temp_sum = temp_sum + expr(i_tissue, i_gene)**2
            end do

            std_dev = sqrt(temp_sum/real(n_replicates, kind=real64))

            if (.not. is_close(std_dev, 0.0_real64)) then
                do concurrent (i_tissue = 1:n_replicates) shared(expr, i_gene, std_dev)
                    expr(i_tissue, i_gene) = expr(i_tissue, i_gene) / std_dev
                end do
            end if
        end do
    end subroutine root_mean_sq_normalization_inplace_helper

    !> AUTHOR_VIVIAN_BASS
    !| Quantile normalization of a gene expression matrix (F42-compliant).
    !| Computes average expression per rank across tissues.
    pure subroutine quantile_normalization(n_genes, n_replicates, expr, normalized_expr, rank_means, tmp_genes_row, tmp_perm, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
        real(real64), dimension(n_replicates, n_genes), intent(out) :: normalized_expr
            !! Normalized `expr`
        real(real64), dimension(n_genes), intent(out) :: rank_means
            !! Preallocated vector to store rank means
        real(real64), dimension(n_genes), intent(out) :: tmp_genes_row
            !! Temporary vector for sorting a tissue in `expr` across genes
        integer(int32), dimension(n_genes), intent(out) :: tmp_perm
            !! Permutation vector
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Error handling
        call set_ok(ierr)

        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_replicates, ierr)

        if (is_err(ierr)) return

        normalized_expr = expr
        call quantile_normalization_inplace_helper(n_genes, n_replicates, normalized_expr, rank_means, tmp_genes_row, tmp_perm)
    end subroutine quantile_normalization

    !> AUTHOR_VIVIAN_BASS
    !| Quantile normalization of a gene expression matrix (F42-compliant).
    !| Computes average expression per rank across tissues.
    pure subroutine quantile_normalization_inplace_helper(n_genes, n_replicates, expr, rank_means, tmp_genes_row, tmp_perm)
        use f42_utils, only: sort_array_heapsort

        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(inout) :: expr
            !! Gene Expression matrix
        real(real64), dimension(n_genes), intent(out) :: rank_means
            !! Preallocated vector to store rank means
        real(real64), dimension(n_genes), intent(out) :: tmp_genes_row
            !! Temporary vector for sorting a tissue in `expr` across genes
        integer(int32), dimension(n_genes), intent(out) :: tmp_perm
            !! Permutation vector

        ! Locals
        integer(int32) :: i_gene, i_tissue

        ! Initialize rank means
        rank_means = 0.0_real64

        do concurrent (i_gene = 1:n_genes) shared(tmp_perm)
            tmp_perm(i_gene) = i_gene
        end do

        ! === First pass: accumulate values by rank across tissues ===
        do i_tissue = 1, n_replicates
            ! Prepare current column and initialize permutation
            do concurrent (i_gene = 1:n_genes) shared(tmp_genes_row, i_tissue, tmp_perm)
                tmp_genes_row(i_gene) = expr(i_tissue, i_gene)
            end do

            ! Sort current column with index tracking
            call sort_array_heapsort(tmp_genes_row, tmp_perm)

            ! Accumulate values for each rank
            do i_gene = 1, n_genes
                rank_means(i_gene) = rank_means(i_gene) + tmp_genes_row(tmp_perm(i_gene))
            end do
        end do

        ! Average the rank values
        do concurrent (i_gene = 1:n_genes) shared(rank_means, n_replicates)
            rank_means(i_gene) = rank_means(i_gene) / real(n_replicates, real64)
        end do

        ! === Second pass: assign averaged values by rank ===
        do i_tissue = 1, n_replicates
            ! Prepare column and reset tmp_permutation
            do concurrent (i_gene = 1:n_genes) shared(tmp_genes_row, i_tissue, tmp_perm)
                tmp_genes_row(i_gene) = expr(i_tissue, i_gene)
            end do

            call sort_array_heapsort(tmp_genes_row, tmp_perm)

            do concurrent (i_gene = 1:n_genes) shared(expr, tmp_perm, rank_means, i_tissue)
                expr(i_tissue, tmp_perm(i_gene)) = rank_means(i_gene)
            end do
        end do
    end subroutine quantile_normalization_inplace_helper

    !> AUTHOR_VIVIAN_BASS
    !| Apply `log2(x + 1)` transformation to each element of the input matrix.
    !| This subroutine performs element-wise `log2(x + 1)` transformation on a
    !| matrix flattened in column-major order. The `log2` is computed via:
    !| `log(x + 1) / log(2)`, which is numerically equivalent and avoids the
    !| non-portable `log2` intrinsic for compatibility with WebAssembly (WASM).
    pure subroutine log2_transformation(n_genes, n_tissues, expr, transformed_expr, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        real(real64), dimension(n_tissues, n_genes), intent(in) :: expr
            !! Gene Expression matrix, from [[tox_normalization(module):calc_tiss_avg(subroutine)]]
        real(real64), dimension(n_tissues, n_genes), intent(out) :: transformed_expr
            !! Log-transformed `expr`
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Error handling
        call set_ok(ierr)

        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_tissues, ierr)

        if (is_err(ierr)) return

        transformed_expr = expr
        call log2_transformation_inplace_helper(n_genes, n_tissues, transformed_expr, ierr)
    end subroutine log2_transformation

    !> AUTHOR_VIVIAN_BASS
    !| Apply `log2(x + 1)` transformation to each element of the input matrix.
    !| This subroutine performs element-wise `log2(x + 1)` transformation on a
    !| matrix flattened in column-major order. The `log2` is computed via:
    !| `log(x + 1) / log(2)`, which is numerically equivalent and avoids the
    !| non-portable `log2` intrinsic for compatibility with WebAssembly (WASM).
    pure subroutine log2_transformation_inplace_helper(n_genes, n_tissues, expr, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        real(real64), dimension(n_tissues, n_genes), intent(inout) :: expr
            !! Gene Expression matrix, from [[tox_normalization(module):calc_tiss_avg(subroutine)]]
        integer(int32), intent(out) :: ierr
            !! Error code
        ! Locals
        integer(int32) :: i_gene, i_group, tmp_ierr
        real(real64) :: expr_val

        call set_ok(ierr)

        ! Loop through all elements in the flattened input matrix
        do concurrent (i_gene = 1:n_genes) shared(n_tissues, expr, ierr)
            do concurrent (i_group = 1:n_tissues) local(tmp_ierr, expr_val) shared(expr, ierr, i_gene)
                ! Apply the log2(x + 1) transformation
                expr_val = expr(i_group, i_gene) + 1.0_real64
                call logx(expr_val, 2.0_real64, expr(i_group, i_gene), tmp_ierr)
                if (is_err(tmp_ierr)) ierr = tmp_ierr
            end do
        end do
    end subroutine log2_transformation_inplace_helper

    !> AUTHOR_VIVIAN_BASS
    !| Calculate tissue averages by averaging replicates within each tissue.
    !| For each tissue of tissue replicates, this subroutine computes the average
    !| expression per gene.
    pure subroutine calc_tiss_avg(n_genes, n_tissues, reps_per_tissue, expr, tissue_averages, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        integer(int32), dimension(n_tissues), intent(in) :: reps_per_tissue
            !! Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
            !! e.g. `[2,3]` means `5` total replicates per gene, with the `expr(1:2, i_gene)` related to the first tissue and `expr(3:, i_gene)` related to the second one.
        real(real64), dimension(sum(reps_per_tissue), n_genes), intent(in) :: expr
            !! Gene Expression matrix
        real(real64), dimension(n_tissues, n_genes), intent(out) :: tissue_averages
            !! Tissue averages per gene
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Error handling
        call set_ok(ierr)

        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_tissues, ierr)
        call validate_all_in_range_int(reps_per_tissue, n_tissues, ierr, min=1_int32)

        if (is_err(ierr)) return

        call calc_tiss_avg_helper(n_genes, n_tissues, reps_per_tissue, expr, tissue_averages)
    end subroutine calc_tiss_avg

    !> AUTHOR_VIVIAN_BASS
    !| (no input validation) Calculate tissue averages by averaging replicates within each group.
    !| For each group of tissue replicates, this subroutine computes the average
    !| expression per gene.
    pure subroutine calc_tiss_avg_helper(n_genes, n_tissues, reps_per_tissue, expr, tissue_averages)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        integer(int32), dimension(n_tissues), intent(in) :: reps_per_tissue
            !! Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
            !! e.g. `[2,3]` means `5` total replicates per gene, with the `expr(1:2, i_gene)` related to the first tissue and `expr(3:, i_gene)` related to the second one.
        real(real64), dimension(sum(reps_per_tissue), n_genes), intent(in) :: expr
            !! Gene Expression matrix
        real(real64), dimension(n_tissues, n_genes), intent(out) :: tissue_averages
            !! Tissue averages per gene

        ! === Local variables ===
        integer(int32) :: i_gene, i_group, i_tissue
        real(real64) :: sum_val
        integer(int32) :: start_idx, stop_idx

        ! === Loop over each group ===
        do concurrent (i_gene = 1:n_genes) shared(n_tissues, reps_per_tissue, expr, tissue_averages)
            do concurrent (i_group = 1:n_tissues) local(start_idx, stop_idx, sum_val) shared(reps_per_tissue, expr, tissue_averages)
                if (i_group == 1) then
                    start_idx = 1
                else
                    start_idx = sum(reps_per_tissue(:i_group-1)) + 1
                end if
                stop_idx = start_idx + reps_per_tissue(i_group) - 1

                sum_val = 0.0_real64
                do concurrent (i_tissue = start_idx:stop_idx) shared(expr, i_gene) reduce(+:sum_val)
                    sum_val = sum_val + expr(i_tissue, i_gene)
                end do

                tissue_averages(i_group, i_gene) = sum_val / real(reps_per_tissue(i_group), real64)
            end do
        end do
    end subroutine calc_tiss_avg_helper

    !> AUTHOR_VIVIAN_BASS
    !| Calculate `log2 fold changes` between condition and control groups.
    !| For each control-condition pair, this subroutine computes the `log2 fold change`
    !| by subtracting the expression value in the control group from the corresponding
    !| value in the condition group, for all genes.
    pure subroutine calc_fchange(n_genes, n_tissues, n_pairs, control_tissues, condition_tissues, expr, fold_changes, ierr)
        ! === Arguments ===
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        integer(int32), intent(in) :: n_pairs
            !! Number of control-condition pairs
        integer(int32), dimension(n_pairs), intent(in) :: control_tissues
            !! Control tissue indices
        integer(int32), dimension(n_pairs), intent(in) :: condition_tissues
            !! Condition tissue indices
        real(real64), dimension(n_tissues, n_genes), intent(in) :: expr
            !! Gene Expression matrix, from [[tox_normalization(module):calc_tiss_avg(subroutine)]]
        real(real64), dimension(n_pairs, n_genes), intent(out) :: fold_changes
            !! Output matrix for fold changes
        integer(int32), intent(out) :: ierr
            !! Error code

        ! === Locals ===
        integer(int32) :: i_gene, i_pair
        integer(int32) :: control_group, cond_group

        ! Error handling
        call set_ok(ierr)

        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_tissues, ierr)
        call validate_dimension_size(n_pairs, ierr)

        if (is_err(ierr)) return

        ! === Loop over each pair ===
        do concurrent (i_gene = 1:n_genes) shared(n_pairs, control_tissues, condition_tissues, expr, fold_changes)
            do concurrent (i_pair = 1:n_pairs) local(control_group, cond_group) shared(i_gene, control_tissues, condition_tissues, expr, fold_changes)
                control_group = control_tissues(i_pair)
                cond_group = condition_tissues(i_pair)
                fold_changes(i_pair, i_gene) = expr(cond_group, i_gene) - expr(control_group, i_gene)
            end do
        end do
    end subroutine calc_fchange

end module tox_normalization

!> C wrapper for [[tox_normalization(module):root_mean_sq_normalization(subroutine)]]
pure subroutine root_mean_sq_normalization_c(n_genes, n_replicates, expr, normalized_expr, ierr) bind(C, name="root_mean_sq_normalization_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_f_pointer, c_loc
    use tox_normalization, only: root_mean_sq_normalization
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes (rows)
    integer(c_int), intent(in), target :: n_replicates
        !! Number of tissues (columns)
    real(c_double), dimension(n_replicates, n_genes), intent(in), target :: expr
        !! Gene Expression matrix
    real(c_double), dimension(n_replicates, n_genes), intent(out), target :: normalized_expr
        !! Normalized `expr`
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_replicates)
    M_CHECK_NON_NULL(expr)
    M_CHECK_NON_NULL(normalized_expr)

    call root_mean_sq_normalization(n_genes, n_replicates, expr, normalized_expr, ierr)

end subroutine root_mean_sq_normalization_c

!> C wrapper for [[tox_normalization(module):normalize_by_std_dev_alloc(subroutine)]]
subroutine normalize_by_std_dev_c(n_genes, n_replicates, expr, normalized_expr, span, degree, ierr) bind(C, name="normalize_by_std_dev_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_normalization, only: normalize_by_std_dev_alloc
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes (rows)
    integer(c_int), intent(in), target :: n_replicates
        !! Number of tissues (columns)
    real(c_double), dimension(n_replicates, n_genes), intent(in), target :: expr
        !! Gene Expression matrix
    real(c_double), dimension(n_replicates, n_genes), intent(out), target :: normalized_expr
        !! Normalized `expr`
    real(c_double), intent(in), target :: span
        !! LOESS span parameter
    integer(c_int), intent(in), target :: degree
        !! LOESS degree parameter
    integer(c_int), intent(out), target:: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_replicates)
    M_CHECK_NON_NULL(expr)
    M_CHECK_NON_NULL(normalized_expr)
    M_CHECK_NON_NULL(span)
    M_CHECK_NON_NULL(degree)

    ! Call the internal routine with 2D arrays
    call normalize_by_std_dev_alloc(n_genes, n_replicates, expr, normalized_expr, span, degree, ierr)

end subroutine normalize_by_std_dev_c

!> C wrapper for [[tox_normalization(module):quantile_normalization(subroutine)]]
pure subroutine quantile_normalization_c(n_genes, n_replicates, expr, normalized_expr, rank_means, tmp_genes_row, tmp_perm, ierr) &
    bind(C, name="quantile_normalization_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_normalization, only: quantile_normalization
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes (rows)
    integer(c_int), intent(in), target :: n_replicates
        !! Number of tissues (columns)
    real(c_double), dimension(n_replicates, n_genes), intent(in), target :: expr
        !! Gene Expression matrix
    real(c_double), dimension(n_replicates, n_genes), intent(out), target :: normalized_expr
        !! Normalized `expr`
    real(c_double), dimension(n_genes), intent(out), target :: tmp_genes_row
        !! Temporary vector for sorting a tissue in `expr` across genes
    real(c_double), dimension(n_genes), intent(out), target :: rank_means
        !! Preallocated vector to store rank means
    integer(c_int), dimension(n_genes), intent(out), target :: tmp_perm
        !! Permutation vector
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_replicates)
    M_CHECK_NON_NULL(expr)
    M_CHECK_NON_NULL(normalized_expr)
    M_CHECK_NON_NULL(tmp_genes_row)
    M_CHECK_NON_NULL(rank_means)
    M_CHECK_NON_NULL(tmp_perm)

    call quantile_normalization(n_genes, n_replicates, expr, normalized_expr, rank_means, tmp_genes_row, tmp_perm, ierr)
end subroutine quantile_normalization_c

!> C wrapper for [[tox_normalization(module):log2_transformation(subroutine)]]
pure subroutine log2_transformation_c(n_genes, n_replicates, expr, transformed_expr, ierr) bind(C, name="log2_transformation_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_normalization, only: log2_transformation
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes (rows)
    integer(c_int), intent(in), target :: n_replicates
        !! Number of tissues (columns)
    real(c_double), dimension(n_replicates, n_genes), intent(in), target :: expr
        !! Gene Expression matrix, from [[tox_normalization(module):calc_tiss_avg(subroutine)]]
    real(c_double), dimension(n_replicates, n_genes), intent(out), target :: transformed_expr
        !! Log-transformed `expr`
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_replicates)
    M_CHECK_NON_NULL(expr)
    M_CHECK_NON_NULL(transformed_expr)

    call log2_transformation(n_genes, n_replicates, expr, transformed_expr, ierr)
end subroutine log2_transformation_c

!> C wrapper for [[tox_normalization(module):calc_tiss_avg(subroutine)]]
pure subroutine calc_tiss_avg_c(n_genes, n_tissues, reps_per_tissue, expr, tissue_averages, ierr) bind(C, name="calc_tiss_avg_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_normalization, only: calc_tiss_avg
    M_USE_NULL_VALIDATION
    implicit none
    integer(c_int), intent(in), target :: n_genes
        !! Number of genes (rows)
    integer(c_int), intent(in), target :: n_tissues
        !! Number of tissues
    integer(c_int), dimension(n_tissues), intent(in), target :: reps_per_tissue
        !! Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
        !! e.g. `[2,3]` means `5` total replicates per gene, with the `expr(1:2, i_gene)` related to the first tissue and `expr(3:, i_gene)` related to the second one.
    real(c_double), dimension(sum(reps_per_tissue), n_genes), intent(in), target :: expr
        !! Gene Expression matrix
    real(c_double), dimension(n_tissues, n_genes), intent(out), target :: tissue_averages
        !! Tissue averages per gene
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_tissues)
    M_CHECK_NON_NULL(reps_per_tissue)
    M_CHECK_NON_NULL(expr)
    M_CHECK_NON_NULL(tissue_averages)

    call calc_tiss_avg(n_genes, n_tissues, reps_per_tissue, expr, tissue_averages, ierr)
end subroutine calc_tiss_avg_c

!> C wrapper for [[tox_normalization(module):calc_fchange(subroutine)]]
pure subroutine calc_fchange_c(n_genes, n_tissues, n_pairs, control_tissues, condition_tissues, expr, fold_changes, ierr) bind(C, name="calc_fchange_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_normalization, only: calc_fchange
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes (rows)
    integer(c_int), intent(in), target :: n_tissues
        !! Number of tissues
    integer(c_int), intent(in), target :: n_pairs
        !! Number of control-condition pairs
    integer(c_int), dimension(n_pairs), intent(in), target :: control_tissues
        !! Control tissue indices
    integer(c_int), dimension(n_pairs), intent(in), target :: condition_tissues
        !! Condition tissue indices
    real(c_double), dimension(n_tissues, n_genes), intent(in), target :: expr
        !! Gene Expression matrix, from [[tox_normalization(module):calc_tiss_avg(subroutine)]]
    real(c_double), dimension(n_pairs, n_genes), intent(out), target :: fold_changes
        !! Output matrix for fold changes
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_tissues)
    M_CHECK_NON_NULL(n_pairs)
    M_CHECK_NON_NULL(control_tissues)
    M_CHECK_NON_NULL(condition_tissues)
    M_CHECK_NON_NULL(expr)
    M_CHECK_NON_NULL(fold_changes)

    call calc_fchange(n_genes, n_tissues, n_pairs, control_tissues, condition_tissues, expr, fold_changes, ierr)
end subroutine calc_fchange_c

!> C wrapper for [[tox_normalization(module):normalization_pipeline_alloc(subroutine)]]
subroutine normalization_pipeline_c(n_genes, n_replicates, expr, log_transformed_expr, reps_per_tissue, n_tissues, span, degree, use_quantile, ierr) bind(C, name="normalization_pipeline_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_normalization, only: normalization_pipeline_alloc
    use tox_conversions, only: c_int_as_logical
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes (rows)
    integer(c_int), intent(in), target :: n_replicates
        !! Number of tissues (columns)
    integer(c_int), intent(in), target :: n_tissues
        !! Number of tissues
    real(c_double), dimension(n_replicates, n_genes), intent(in), target :: expr
        !! Gene Expression matrix
    integer(c_int), dimension(n_tissues), intent(in), target :: reps_per_tissue
        !! Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
        !! e.g. `[2,3]` means `5` total replicates per gene, with the `expr(1:2, i_gene)` related to the first tissue and `expr(3:, i_gene)` related to the second one.
    real(c_double), dimension(n_tissues, n_genes), intent(out), target :: log_transformed_expr
        !! Log-transformed grouped `expr`
    real(c_double), intent(in), target :: span
        !! LOESS span parameter
    integer(c_int), intent(in), target :: degree
        !! LOESS degree parameter
    integer(c_int), intent(in), target :: use_quantile
        !! Use quantile normalization (0/1 logical)
    integer(c_int), intent(out), target :: ierr
        !! Error code

    logical :: use_quantile_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_replicates)
    M_CHECK_NON_NULL(n_tissues)
    M_CHECK_NON_NULL(expr)
    M_CHECK_NON_NULL(reps_per_tissue)
    M_CHECK_NON_NULL(log_transformed_expr)
    M_CHECK_NON_NULL(span)
    M_CHECK_NON_NULL(degree)
    M_CHECK_NON_NULL(use_quantile)

    call c_int_as_logical(use_quantile, use_quantile_f)

    call normalization_pipeline_alloc(n_genes, n_replicates, expr, log_transformed_expr, reps_per_tissue, n_tissues, span, degree, use_quantile_f, ierr)
end subroutine normalization_pipeline_c

!> C wrapper for [[tox_normalization(module):normalize_unit_length(subroutine)]]
pure subroutine normalize_unit_length_c(vector, n_dims, ierr) bind(C, name="normalize_unit_length_c")
    use tox_normalization, only: normalize_unit_length
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_dims
        !! number of elements in `vector`
    real(c_double), dimension(n_dims), intent(inout), target :: vector
        !! Vector that will be normalized to unit length
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_dims)
    M_CHECK_NON_NULL(vector)

    call normalize_unit_length(vector, n_dims, ierr)
end subroutine normalize_unit_length_c
