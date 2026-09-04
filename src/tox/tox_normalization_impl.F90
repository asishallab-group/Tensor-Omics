#include <src/macros.h>

!> Normalization of expression data: putting samples on a comparable scale before anything
!| distance-based is computed from them.
!|
!| Several schemes, each its own entry point rather than a mode: a log2 transformation, scaling
!| by standard deviation (LOESS-smoothed against expression level, so the correction follows the
!| mean-variance trend rather than assuming one), root-mean-square scaling, quantile
!| normalization, and normalization to unit length.
!|
!| `normalization_pipeline` chains them in the order the analysis expects, and is what most
!| callers want. The individual routines are published for a caller assembling their own.
module tox_normalization_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err, ERR_DIVISION_BY_ZERO, ERR_INVALID_INPUT, is_err, &
                          validate_in_range_real, validate_all_in_range_real, ERR_SIZE_MISMATCH
    use f42_math_impl, only: is_close, logx_helper, above, mean, std_dev
    use f42_vector_impl, only: norm
    use tox_loess_impl, only: loess_fit_robust_impl

#define CM_LOESS_SPAN_DEFAULT 0.7_real64
#define CM_LOESS_DEGREE_DEFAULT 2_int32
#define CM_LOESS_ROBUST_ITERS 3_int32

contains

    !> summary: Normalizes an input vector to unit length in-place
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine normalize_unit_length_impl(vector, n_dims, ierr)
        integer(int32), intent(in) :: n_dims
            !! number of elements in `vector`
        real(real64), dimension(n_dims), intent(inout) :: vector
            !! Vector that will be normalized to unit length
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_dim
        real(real64) :: vector_norm

        ! The norm is a derived quantity, so its zero / non-finite guards are runtime checks here.
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
    end subroutine normalize_unit_length_impl

    !> summary: Complete normalization pipeline for gene expression data.
    !| AUTHOR_VIVIAN_BASS
    !| Performs: std dev normalization, quantile normalization, replicate averaging, log2(x+1) transformation.
    !| Final result is in log_transformed_expr. If fold change is needed, call calc_fchange separately.
    subroutine normalization_pipeline_impl(n_genes, n_replicates, expr, log_transformed_expr, reps_per_tissue, n_tissues, &
                                             tmp_expr_copy, tmp_loess_y, tmp_indices_used, tmp_yhat_global, &
                                             tmp_int_workspace, int_workspace_size, tmp_real_workspace, real_workspace_size, &
                                             tmp_hat_diag, tmp_loess_weights, tmp_eval_points, &
                                             tmp_robust_weights, tmp_combined_weights, tmp_residuals, tmp_permutation_indices, &
                                             span, degree, use_quantile, ierr)

        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        integer(int32), dimension(n_tissues), intent(in) :: reps_per_tissue
            !! Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
            !! e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
        real(real64), dimension(n_tissues, n_genes), intent(out), target :: log_transformed_expr
            !! Log-transformed grouped `expr`

        ! The pipeline runs in place, so it works on a copy: `expr` is the caller's.
        real(real64), dimension(n_replicates, n_genes), intent(out) :: tmp_expr_copy
            !! Work matrix the pipeline normalizes in place
        ! The (mean, sd) scatter the LOESS fit is built from. Only the vectors that cannot be
        ! taken from the output buffer below are passed in; see the aliasing note in the body.
        real(real64), dimension(n_genes), intent(out), target :: tmp_loess_y
            !! Work vector for the empirical standard deviations (Y-axis for LOESS)
        integer(int32), dimension(n_genes), intent(out) :: tmp_indices_used
            !! Work vector mapping the fitted points back to gene indices
        real(real64), dimension(n_genes), intent(out), target :: tmp_yhat_global
            !! Work vector for the fitted standard deviations (LOESS predictions)

        ! LOESS workspace. Sized for `n_genes`, the most points the fit can ever see: the
        ! genes that carry no variance are dropped from it, never added to it.
        integer(int32), intent(in) :: int_workspace_size
            !! Length of integer workspace.
            !! DM_OUTPUT_FROM(int_workspace_size, tox_loess_required_workspace, tox_loess_impl, AUTO)
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_genes     |
            !! | save_factorization    | .false.     |
        integer(int32), dimension(int_workspace_size), intent(out) :: tmp_int_workspace
            !! Integer workspace array
        integer(int32), intent(in) :: real_workspace_size
            !! Length of real workspace.
            !! DM_OUTPUT_FROM(real_workspace_size, tox_loess_required_workspace, tox_loess_impl, AUTO)
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_genes     |
            !! | save_factorization    | .false.     |
        real(real64), dimension(real_workspace_size), intent(out) :: tmp_real_workspace
            !! Real workspace array
        real(real64), dimension(n_genes), intent(out) :: tmp_hat_diag
            !! Diagonal elements of the LOESS hat matrix
        real(real64), dimension(n_genes), intent(out) :: tmp_loess_weights
            !! Per-point weights handed to the LOESS fit
        real(real64), dimension(n_genes, 1), intent(out) :: tmp_eval_points
            !! Points the fitted curve is evaluated at
        real(real64), dimension(n_genes), intent(out) :: tmp_robust_weights
            !! Robust bisquare weights of the LOESS fit
        real(real64), dimension(n_genes), intent(out) :: tmp_combined_weights
            !! Combined weights of the LOESS fit
        real(real64), dimension(n_genes), intent(out) :: tmp_residuals
            !! Residuals of the LOESS fit
        integer(int32), dimension(n_genes), intent(out) :: tmp_permutation_indices
            !! Permutation indices of the LOESS fit

        real(real64), intent(in), optional :: span
            !! LOESS span parameter.
            !! DM_DEFAULT(CM_LOESS_SPAN_DEFAULT)
        integer(int32), intent(in), optional :: degree
            !! LOESS degree parameter.
            !! DM_DEFAULT(CM_LOESS_DEGREE_DEFAULT)
        logical(c_bool), intent(in), optional :: use_quantile
            !! Use quantile normalization.
            !! DM_DEFAULT(.false.)
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Local variables
        logical(c_bool) :: actual_use_quantile

        real(real64), dimension(:, :), pointer :: log_transformed_expr_transposed_view
        real(real64), dimension(:), pointer :: tmp_loess_x_ptr, tmp_loess_y_ptr, tmp_yhat_global_ptr

        ! Error handling
        call set_ok(ierr)

        ! sum(reps_per_tissue) must equal n_replicates -- a relation between arguments.
        if (sum(reps_per_tissue) /= n_replicates) then
            call set_err(ierr, ERR_SIZE_MISMATCH)
            return
        end if

        M_DEFAULT_VAL(use_quantile, actual_use_quantile, .false.)

        ! Reuse spare columns of the (n_genes, n_tissues) output buffer as scratch space for the
        ! per-gene LOESS x/y/yhat vectors (each length n_genes) below, instead of reading the
        ! passed-in work vectors, since those columns are still unwritten at this point and get
        ! fully overwritten by calc_tiss_avg_helper further down before they are read as output.
        ! Only the columns that exist (n_tissues >= 2 for column 2, >= 3 for column 3) can be
        ! reused this way; the remaining vectors come from the caller. Column 1 always exists, so
        ! the x vector is never a dummy at all.
        log_transformed_expr_transposed_view(1:n_genes, 1:n_tissues) => log_transformed_expr
        tmp_loess_x_ptr => log_transformed_expr_transposed_view(:, 1)
        select case (n_tissues)
            case (1)
                tmp_loess_y_ptr => tmp_loess_y
                tmp_yhat_global_ptr => tmp_yhat_global
            case (2)
                tmp_loess_y_ptr => log_transformed_expr_transposed_view(:, 2)
                tmp_yhat_global_ptr => tmp_yhat_global
            case (3:)
                tmp_loess_y_ptr => log_transformed_expr_transposed_view(:, 2)
                tmp_yhat_global_ptr => log_transformed_expr_transposed_view(:, 3)
        end select

        tmp_expr_copy = expr

        ! Step 1: Normalize per-gene by std dev
        call normalize_by_std_dev_inplace_helper(n_genes, n_replicates, tmp_expr_copy, &
                                    tmp_loess_x_ptr, tmp_loess_y_ptr, tmp_indices_used, tmp_yhat_global_ptr, &
                                    tmp_int_workspace, int_workspace_size, tmp_real_workspace, real_workspace_size, &
                                    tmp_hat_diag, tmp_loess_weights, tmp_eval_points, &
                                    tmp_robust_weights, tmp_combined_weights, tmp_residuals, tmp_permutation_indices, &
                                    span, degree, ierr)
        if (is_err(ierr)) return

        ! Step 2: Quantile normalization (conditional). The three LOESS vectors are spent by now,
        ! so they serve as the rank/row/permutation scratch here.
        if (actual_use_quantile) then
            call quantile_normalization_inplace_helper(n_genes, n_replicates, tmp_expr_copy, &
                                    tmp_loess_x_ptr, tmp_loess_y_ptr, tmp_indices_used)
        end if

        call calc_tiss_avg_helper(n_genes, n_tissues, reps_per_tissue, tmp_expr_copy, log_transformed_expr)

        ! Step 4: Log2(x+1) transformation
        call log2_transformation_inplace_helper(n_genes, n_tissues, log_transformed_expr, ierr)
        if (is_err(ierr)) return
    end subroutine normalization_pipeline_impl

    !> summary: Normalizes each gene's expression vector using LOESS-stabilized standard deviation.
    !| AUTHOR_VIVIAN_BASS
    !| This procedure applies a global stabilization based on the relationship between
    !| gene-wise mean expression and empirical standard deviation.
    subroutine normalize_by_std_dev_impl(n_genes, n_replicates, expr, normalized_expr, &
                                    tmp_loess_x, tmp_loess_y, tmp_indices_used, tmp_yhat_global, &
                                    tmp_int_workspace, int_workspace_size, tmp_real_workspace, real_workspace_size, &
                                    tmp_hat_diag, tmp_loess_weights, tmp_eval_points, &
                                    tmp_robust_weights, tmp_combined_weights, tmp_residuals, tmp_permutation_indices, &
                                    span, degree, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        real(real64), dimension(n_replicates, n_genes), intent(out) :: normalized_expr
            !! Normalized `expr`

        ! The (mean, sd) scatter the LOESS fit is built from
        real(real64), dimension(n_genes), intent(out) :: tmp_loess_x
            !! Work vector for the mean values (X-axis for LOESS)
        real(real64), dimension(n_genes), intent(out) :: tmp_loess_y
            !! Work vector for the empirical standard deviations (Y-axis for LOESS)
        integer(int32), dimension(n_genes), intent(out) :: tmp_indices_used
            !! Work vector mapping the fitted points back to gene indices
        real(real64), dimension(n_genes), intent(out) :: tmp_yhat_global
            !! Work vector for the fitted standard deviations (LOESS predictions)

        ! LOESS workspace. Sized for `n_genes`, the most points the fit can ever see: the
        ! genes that carry no variance are dropped from it, never added to it.
        integer(int32), intent(in) :: int_workspace_size
            !! Length of integer workspace.
            !! DM_OUTPUT_FROM(int_workspace_size, tox_loess_required_workspace, tox_loess_impl, AUTO)
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_genes     |
            !! | save_factorization    | .false.     |
        integer(int32), dimension(int_workspace_size), intent(out) :: tmp_int_workspace
            !! Integer workspace array
        integer(int32), intent(in) :: real_workspace_size
            !! Length of real workspace.
            !! DM_OUTPUT_FROM(real_workspace_size, tox_loess_required_workspace, tox_loess_impl, AUTO)
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_genes     |
            !! | save_factorization    | .false.     |
        real(real64), dimension(real_workspace_size), intent(out) :: tmp_real_workspace
            !! Real workspace array
        real(real64), dimension(n_genes), intent(out) :: tmp_hat_diag
            !! Diagonal elements of the LOESS hat matrix
        real(real64), dimension(n_genes), intent(out) :: tmp_loess_weights
            !! Per-point weights handed to the LOESS fit
        real(real64), dimension(n_genes, 1), intent(out) :: tmp_eval_points
            !! Points the fitted curve is evaluated at
        real(real64), dimension(n_genes), intent(out) :: tmp_robust_weights
            !! Robust bisquare weights of the LOESS fit
        real(real64), dimension(n_genes), intent(out) :: tmp_combined_weights
            !! Combined weights of the LOESS fit
        real(real64), dimension(n_genes), intent(out) :: tmp_residuals
            !! Residuals of the LOESS fit
        integer(int32), dimension(n_genes), intent(out) :: tmp_permutation_indices
            !! Permutation indices of the LOESS fit

        real(real64), intent(in), optional :: span
            !! LOESS span parameter.
            !! DM_DEFAULT(CM_LOESS_SPAN_DEFAULT)
        integer(int32), intent(in), optional :: degree
            !! LOESS degree parameter.
            !! DM_DEFAULT(CM_LOESS_DEGREE_DEFAULT)
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        normalized_expr = expr
        call normalize_by_std_dev_inplace_helper(n_genes, n_replicates, normalized_expr, &
                                    tmp_loess_x, tmp_loess_y, tmp_indices_used, tmp_yhat_global, &
                                    tmp_int_workspace, int_workspace_size, tmp_real_workspace, real_workspace_size, &
                                    tmp_hat_diag, tmp_loess_weights, tmp_eval_points, &
                                    tmp_robust_weights, tmp_combined_weights, tmp_residuals, tmp_permutation_indices, &
                                    span, degree, ierr)
    end subroutine normalize_by_std_dev_impl

    !> AUTHOR_VIVIAN_BASS
    !| Normalizes each gene's expression vector using LOESS-stabilized standard deviation.
    !| This procedure applies a global stabilization based on the relationship between
    !| gene-wise mean expression and empirical standard deviation.
    subroutine normalize_by_std_dev_inplace_helper(n_genes, n_replicates, expr, &
                                    tmp_loess_x, tmp_loess_y, tmp_indices_used, tmp_yhat_global, &
                                    tmp_int_workspace, int_workspace_size, tmp_real_workspace, real_workspace_size, &
                                    tmp_hat_diag, tmp_loess_weights, tmp_eval_points, &
                                    tmp_robust_weights, tmp_combined_weights, tmp_residuals, tmp_permutation_indices, &
                                    span, degree, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(inout) :: expr
            !! Gene Expression matrix

        ! Buffers for LOESS fitting (owned by the caller, so this routine allocates nothing)
        real(real64), dimension(n_genes), intent(out) :: tmp_loess_x
            !! Mean values (X-axis for LOESS)
        real(real64), dimension(n_genes), intent(out) :: tmp_loess_y
            !! Empirical standard deviation values (Y-axis for LOESS)
        integer(int32), dimension(n_genes), intent(out) :: tmp_indices_used
            !! Mapping back to gene indices
        real(real64), dimension(n_genes), intent(out) :: tmp_yhat_global
            !! Fitted standard deviation values (LOESS predictions)

        integer(int32), intent(in) :: int_workspace_size
            !! Length of integer workspace
        integer(int32), dimension(int_workspace_size), intent(out) :: tmp_int_workspace
            !! Integer workspace array
        integer(int32), intent(in) :: real_workspace_size
            !! Length of real workspace
        real(real64), dimension(real_workspace_size), intent(out) :: tmp_real_workspace
            !! Real workspace array
        real(real64), dimension(n_genes), intent(out) :: tmp_hat_diag
            !! Diagonal elements of the LOESS hat matrix
        real(real64), dimension(n_genes), intent(out) :: tmp_loess_weights
            !! Per-point weights handed to the LOESS fit
        real(real64), dimension(n_genes, 1), intent(out) :: tmp_eval_points
            !! Points the fitted curve is evaluated at
        real(real64), dimension(n_genes), intent(out) :: tmp_robust_weights
            !! Robust bisquare weights of the LOESS fit
        real(real64), dimension(n_genes), intent(out) :: tmp_combined_weights
            !! Combined weights of the LOESS fit
        real(real64), dimension(n_genes), intent(out) :: tmp_residuals
            !! Residuals of the LOESS fit
        integer(int32), dimension(n_genes), intent(out) :: tmp_permutation_indices
            !! Permutation indices of the LOESS fit

        real(real64), intent(in), optional :: span
            !! LOESS span parameter.
            !! DM_DEFAULT(CM_LOESS_SPAN_DEFAULT)
        integer(int32), intent(in), optional :: degree
            !! LOESS degree parameter.
            !! DM_DEFAULT(CM_LOESS_DEGREE_DEFAULT)
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

            ! Genes with zero variance across replicates carry no information about the mean-vs-sd
            ! trend and would only produce a degenerate (division-by-zero) target for LOESS, so they
            ! are dropped from the fit here. Since n_valid <= i_gene always, this compacts the arrays
            ! in place (overwriting already-consumed slots) rather than needing a separate buffer.
            if (is_close(tmp_loess_y(i_gene), 0.0_real64)) cycle

            n_valid = n_valid + 1
            tmp_loess_x(n_valid) = tmp_loess_x(i_gene)
            tmp_loess_y(n_valid) = tmp_loess_y(i_gene)
            tmp_indices_used(n_valid) = i_gene
        end do

        ! Require a handful of points so the LOESS fit below is not driven by noise from too few
        ! (mean, sd) pairs.
        if (n_valid < 5) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        ! Fit the mean-vs-sd trend robustly. This is the expert tier of the LOESS API, so the
        ! workspace, the uniform weights and the evaluation points are prepared here: the fit
        ! is evaluated at its own training points, every point weighted equally, and the
        ! neighbourhood may span the whole sample. The degenerate cases (too few points, a
        ! near-constant x) the fit answers itself.
        tmp_loess_weights(1:n_valid) = 1.0_real64
        tmp_eval_points(1:n_valid, 1) = tmp_loess_x(1:n_valid)

        call loess_fit_robust_impl(n_valid, tmp_loess_x(1:n_valid), tmp_loess_y(1:n_valid), &
                                     tmp_loess_weights(1:n_valid), tmp_eval_points(1:n_valid, 1), &
                                     actual_span, actual_degree, &
                                     n_valid, .false._c_bool, .false._c_bool, CM_LOESS_ROBUST_ITERS, &
                                     tmp_int_workspace, int_workspace_size, &
                                     tmp_real_workspace, real_workspace_size, tmp_hat_diag(1:n_valid), &
                                     tmp_robust_weights(1:n_valid), tmp_combined_weights(1:n_valid), &
                                     tmp_residuals(1:n_valid), tmp_permutation_indices(1:n_valid), &
                                     tmp_yhat_global(1:n_valid), ierr)
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

    !> summary: Normalizes each gene's expression vector using `sqrt(mean(x^2))`
    !| AUTHOR_VIVIAN_BASS
    !| across tissues (not classical standard deviation).
    pure subroutine root_mean_sq_normalization_impl(n_genes, n_replicates, expr, normalized_expr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        real(real64), dimension(n_replicates, n_genes), intent(out) :: normalized_expr
            !! Normalized `expr`

        normalized_expr = expr
        call root_mean_sq_normalization_inplace_helper(n_genes, n_replicates, normalized_expr)
    end subroutine root_mean_sq_normalization_impl

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

    !> summary: Quantile normalization of a gene expression matrix (F42-compliant).
    !| AUTHOR_VIVIAN_BASS
    !| Computes average expression per rank across tissues.
    pure subroutine quantile_normalization_impl(n_genes, n_replicates, expr, normalized_expr, rank_means, tmp_genes_row, tmp_perm)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        real(real64), dimension(n_replicates, n_genes), intent(out) :: normalized_expr
            !! Normalized `expr`
        real(real64), dimension(n_genes), intent(out) :: rank_means
            !! The mean of each rank across tissues, one per gene
        real(real64), dimension(n_genes), intent(out) :: tmp_genes_row
            !! Temporary vector for sorting a tissue in `expr` across genes
        integer(int32), dimension(n_genes), intent(out) :: tmp_perm
            !! Permutation vector

        normalized_expr = expr
        call quantile_normalization_inplace_helper(n_genes, n_replicates, normalized_expr, rank_means, tmp_genes_row, tmp_perm)
    end subroutine quantile_normalization_impl

    !> AUTHOR_VIVIAN_BASS
    !| Quantile normalization of a gene expression matrix (F42-compliant).
    !| Computes average expression per rank across tissues.
    pure subroutine quantile_normalization_inplace_helper(n_genes, n_replicates, expr, rank_means, tmp_genes_row, tmp_perm)
        use f42_sort_impl, only: sort_array_heapsort

        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(inout) :: expr
            !! Gene Expression matrix
        real(real64), dimension(n_genes), intent(out) :: rank_means
            !! The mean of each rank across tissues, one per gene
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

    !> summary: Apply `log2(x + 1)` transformation to each element of the input matrix.
    !| AUTHOR_VIVIAN_BASS
    !| This subroutine performs element-wise `log2(x + 1)` transformation on a
    !| matrix flattened in column-major order. The `log2` is computed via:
    !| `log(x + 1) / log(2)`, which is numerically equivalent and avoids the
    !| non-portable `log2` intrinsic for compatibility with WebAssembly (WASM).
    pure subroutine log2_transformation_impl(n_genes, n_tissues, expr, transformed_expr, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        real(real64), dimension(n_tissues, n_genes), intent(in) :: expr
            !! Gene Expression matrix, from [[tox_normalization(module):calc_tiss_avg(subroutine)]]
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        real(real64), dimension(n_tissues, n_genes), intent(out) :: transformed_expr
            !! Log-transformed `expr`
        integer(int32), intent(out) :: ierr
            !! Error code

        transformed_expr = expr
        call log2_transformation_inplace_helper(n_genes, n_tissues, transformed_expr, ierr)
    end subroutine log2_transformation_impl

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
        integer(int32) :: i_gene, i_group
        real(real64) :: expr_val

        call set_ok(ierr)

        ! Validate every element up front (sequentially, as `ierr` is a shared scalar here) so that
        ! `log2(x + 1)` is only ever evaluated for `x + 1 > 0`, i.e. `x > -1`. With the inputs
        ! guaranteed valid, the transformation itself can run as a race-free `do concurrent` calling
        ! the non-validating `logx_helper` -- no per-iteration write to the shared `ierr` is needed.
        call validate_all_in_range_real(expr, n_genes*n_tissues, ierr, min=above(-1.0_real64))
        if (is_err(ierr)) return

        ! Apply the log2(x + 1) transformation to every element in the flattened input matrix
        do concurrent(i_gene=1:n_genes) shared(expr, n_tissues)
            do concurrent(i_group=1:n_tissues) local(expr_val) shared(expr, i_gene)
                expr_val = expr(i_group, i_gene) + 1.0_real64
                call logx_helper(expr_val, 2.0_real64, expr(i_group, i_gene))
            end do
        end do
    end subroutine log2_transformation_inplace_helper

    !> summary: Calculate tissue averages by averaging replicates within each tissue.
    !| AUTHOR_VIVIAN_BASS
    !| For each tissue of tissue replicates, this subroutine computes the average
    !| expression per gene.
    pure subroutine calc_tiss_avg_impl(n_genes, n_tissues, reps_per_tissue, expr, tissue_averages)
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        integer(int32), dimension(n_tissues), intent(in) :: reps_per_tissue
            !! Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
            !! e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
            !! DM_MIN(1_int32)
        real(real64), dimension(sum(reps_per_tissue), n_genes), intent(in) :: expr
            !! Gene Expression matrix
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        real(real64), dimension(n_tissues, n_genes), intent(out) :: tissue_averages
            !! Tissue averages per gene

        call calc_tiss_avg_helper(n_genes, n_tissues, reps_per_tissue, expr, tissue_averages)
    end subroutine calc_tiss_avg_impl

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
            !! e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
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
            start_idx = 1
            do i_group = 1, n_tissues
                stop_idx = start_idx + reps_per_tissue(i_group) - 1

                sum_val = 0.0_real64
                do concurrent (i_tissue = start_idx:stop_idx) shared(expr, i_gene) reduce(+:sum_val)
                    sum_val = sum_val + expr(i_tissue, i_gene)
                end do

                tissue_averages(i_group, i_gene) = sum_val / real(reps_per_tissue(i_group), real64)
                start_idx = stop_idx + 1
            end do
        end do
    end subroutine calc_tiss_avg_helper

    !> summary: Calculate `log2 fold changes` between condition and control groups.
    !| AUTHOR_VIVIAN_BASS
    !| For each control-condition pair, this subroutine computes the `log2 fold change`
    !| by subtracting the expression value in the control group from the corresponding
    !| value in the condition group, for all genes.
    pure subroutine calc_fchange_impl(n_genes, n_tissues, n_pairs, control_tissues, condition_tissues, expr, fold_changes)
        ! === Arguments ===
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        integer(int32), intent(in) :: n_pairs
            !! Number of control-condition pairs
        integer(int32), dimension(n_pairs), intent(in) :: control_tissues
            !! Control tissue indices
            !! DM_MIN(1_int32)
            !! DM_MAX(n_tissues)
        integer(int32), dimension(n_pairs), intent(in) :: condition_tissues
            !! Condition tissue indices
            !! DM_MIN(1_int32)
            !! DM_MAX(n_tissues)
        real(real64), dimension(n_tissues, n_genes), intent(in) :: expr
            !! Gene Expression matrix, from [[tox_normalization(module):calc_tiss_avg(subroutine)]]
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        real(real64), dimension(n_pairs, n_genes), intent(out) :: fold_changes
            !! Output matrix for fold changes

        ! === Locals ===
        integer(int32) :: i_gene, i_pair
        integer(int32) :: control_group, cond_group

        ! === Loop over each pair ===
        do concurrent (i_gene = 1:n_genes) shared(n_pairs, control_tissues, condition_tissues, expr, fold_changes)
            do concurrent (i_pair = 1:n_pairs) local(control_group, cond_group) shared(i_gene, control_tissues, condition_tissues, expr, fold_changes)
                control_group = control_tissues(i_pair)
                cond_group = condition_tissues(i_pair)
                fold_changes(i_pair, i_gene) = expr(cond_group, i_gene) - expr(control_group, i_gene)
            end do
        end do
    end subroutine calc_fchange_impl

end module tox_normalization_impl
