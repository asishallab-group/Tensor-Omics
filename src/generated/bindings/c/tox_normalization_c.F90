#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_normalization(module)]]
!| Normalization of expression data: putting samples on a comparable scale before anything
!| distance-based is computed from them.
!|
!| Several schemes, each its own entry point rather than a mode: a log2 transformation, scaling
!| by standard deviation (LOESS-smoothed against expression level, so the correction follows the
!| mean-variance trend rather than assuming one), root-mean-square scaling, quantile
!| normalization, and normalization to unit length.
!|
!| `normalization_pipeline` chains them in the order the analysis expects, and is what most
!| callers want. The individual routines are published for a caller assembling their own.
module tox_normalization_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: normalize_unit_length_c
    public :: normalization_pipeline_c
    public :: normalization_pipeline_expert_c
    public :: normalize_by_std_dev_c
    public :: normalize_by_std_dev_expert_c
    public :: root_mean_sq_normalization_c
    public :: quantile_normalization_c
    public :: quantile_normalization_expert_c
    public :: log2_transformation_c
    public :: calc_tiss_avg_c
    public :: calc_fchange_c

contains

    !> summary: C-wrapper for [[tox_normalization(module):normalize_unit_length(subroutine)]]
    subroutine normalize_unit_length_c(&
            vector,&
            n_dims,&
            ierr&
        ) bind(C, name="normalize_unit_length_c")
        use tox_normalization, only: normalize_unit_length

        integer(c_int), intent(in), target :: n_dims
            !! number of elements in `vector`
        real(c_double), dimension(n_dims), intent(inout), target :: vector
            !! Vector that will be normalized to unit length
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_ARRAY_NON_NULL(vector, n_dims)

        call normalize_unit_length(&
            vector = vector,&
            n_dims = n_dims,&
            ierr = ierr&
        )
    end subroutine normalize_unit_length_c

    !> summary: C-wrapper for [[tox_normalization(module):normalization_pipeline(subroutine)]]
    !| Final result is in log_transformed_expr. If fold change is needed, call calc_fchange separately.
    subroutine normalization_pipeline_c(&
            n_genes,&
            n_replicates,&
            expr,&
            log_transformed_expr,&
            reps_per_tissue,&
            n_tissues,&
            span,&
            degree,&
            use_quantile,&
            ierr&
        ) bind(C, name="normalization_pipeline_c")
        use tox_normalization, only: normalization_pipeline

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes (rows)
        integer(c_int), intent(in), target :: n_replicates
            !! Number of replicates per gene
        integer(c_int), intent(in), target :: n_tissues
            !! Number of tissues
        real(c_double), dimension(n_replicates, n_genes), intent(in), target :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_tissues, n_genes), intent(out), target :: log_transformed_expr
            !! Log-transformed grouped `expr`
        integer(c_int), dimension(n_tissues), intent(in), target :: reps_per_tissue
            !! Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
            !! e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
        real(c_double), intent(in), target :: span
            !! LOESS span parameter.
            !! The default value is `0.7_real64`.
        integer(c_int), intent(in), target :: degree
            !! LOESS degree parameter.
            !! The default value is `2_int32`.
        logical(c_bool), intent(in), target :: use_quantile
            !! Use quantile normalization.
            !! The default value is `.false.`.
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical :: use_quantile_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_replicates)
        M_CHECK_NON_NULL(n_tissues)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_NON_NULL(use_quantile)
        M_CHECK_ARRAY_NON_NULL(expr, n_replicates * n_genes)
        M_CHECK_ARRAY_NON_NULL(log_transformed_expr, n_tissues * n_genes)
        M_CHECK_ARRAY_NON_NULL(reps_per_tissue, n_tissues)

        use_quantile_f = use_quantile

        call normalization_pipeline(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            log_transformed_expr = log_transformed_expr,&
            reps_per_tissue = reps_per_tissue,&
            n_tissues = n_tissues,&
            span = span,&
            degree = degree,&
            use_quantile = use_quantile_f,&
            ierr = ierr&
        )
    end subroutine normalization_pipeline_c

    !> summary: C-wrapper for [[tox_normalization(module):normalization_pipeline_expert(subroutine)]]
    !| Final result is in log_transformed_expr. If fold change is needed, call calc_fchange separately.
    subroutine normalization_pipeline_expert_c(&
            n_genes,&
            n_replicates,&
            expr,&
            log_transformed_expr,&
            reps_per_tissue,&
            n_tissues,&
            tmp_expr_copy,&
            tmp_loess_y,&
            tmp_indices_used,&
            tmp_yhat_global,&
            tmp_int_workspace,&
            int_workspace_size,&
            tmp_real_workspace,&
            real_workspace_size,&
            tmp_hat_diag,&
            tmp_loess_weights,&
            tmp_eval_points,&
            tmp_robust_weights,&
            tmp_combined_weights,&
            tmp_residuals,&
            tmp_permutation_indices,&
            span,&
            degree,&
            use_quantile,&
            ierr&
        ) bind(C, name="normalization_pipeline_expert_c")
        use tox_normalization, only: normalization_pipeline_expert

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes (rows)
        integer(c_int), intent(in), target :: n_replicates
            !! Number of replicates per gene
        integer(c_int), intent(in), target :: n_tissues
            !! Number of tissues
        integer(c_int), intent(in), target :: int_workspace_size
            !! Length of integer workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `int_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_genes     |
            !! | save_factorization    | .false.     |
        integer(c_int), intent(in), target :: real_workspace_size
            !! Length of real workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `real_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_genes     |
            !! | save_factorization    | .false.     |
        real(c_double), dimension(n_replicates, n_genes), intent(in), target :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_tissues, n_genes), intent(out), target :: log_transformed_expr
            !! Log-transformed grouped `expr`
        integer(c_int), dimension(n_tissues), intent(in), target :: reps_per_tissue
            !! Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
            !! e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
        real(c_double), dimension(n_replicates, n_genes), intent(out), target :: tmp_expr_copy
            !! Work matrix the pipeline normalizes in place
        real(c_double), dimension(n_genes), intent(out), target :: tmp_loess_y
            !! Work vector for the empirical standard deviations (Y-axis for LOESS)
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_indices_used
            !! Work vector mapping the fitted points back to gene indices
        real(c_double), dimension(n_genes), intent(out), target :: tmp_yhat_global
            !! Work vector for the fitted standard deviations (LOESS predictions)
        integer(c_int), dimension(int_workspace_size), intent(out), target :: tmp_int_workspace
            !! Integer workspace array
        real(c_double), dimension(real_workspace_size), intent(out), target :: tmp_real_workspace
            !! Real workspace array
        real(c_double), dimension(n_genes), intent(out), target :: tmp_hat_diag
            !! Diagonal elements of the LOESS hat matrix
        real(c_double), dimension(n_genes), intent(out), target :: tmp_loess_weights
            !! Per-point weights handed to the LOESS fit
        real(c_double), dimension(n_genes, 1), intent(out), target :: tmp_eval_points
            !! Points the fitted curve is evaluated at
        real(c_double), dimension(n_genes), intent(out), target :: tmp_robust_weights
            !! Robust bisquare weights of the LOESS fit
        real(c_double), dimension(n_genes), intent(out), target :: tmp_combined_weights
            !! Combined weights of the LOESS fit
        real(c_double), dimension(n_genes), intent(out), target :: tmp_residuals
            !! Residuals of the LOESS fit
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_permutation_indices
            !! Permutation indices of the LOESS fit
        real(c_double), intent(in), target :: span
            !! LOESS span parameter.
            !! The default value is `0.7_real64`.
        integer(c_int), intent(in), target :: degree
            !! LOESS degree parameter.
            !! The default value is `2_int32`.
        logical(c_bool), intent(in), target :: use_quantile
            !! Use quantile normalization.
            !! The default value is `.false.`.
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical :: use_quantile_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_replicates)
        M_CHECK_NON_NULL(n_tissues)
        M_CHECK_NON_NULL(int_workspace_size)
        M_CHECK_NON_NULL(real_workspace_size)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_NON_NULL(use_quantile)
        M_CHECK_ARRAY_NON_NULL(expr, n_replicates * n_genes)
        M_CHECK_ARRAY_NON_NULL(log_transformed_expr, n_tissues * n_genes)
        M_CHECK_ARRAY_NON_NULL(reps_per_tissue, n_tissues)
        M_CHECK_ARRAY_NON_NULL(tmp_expr_copy, n_replicates * n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_loess_y, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_indices_used, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_yhat_global, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_int_workspace, int_workspace_size)
        M_CHECK_ARRAY_NON_NULL(tmp_real_workspace, real_workspace_size)
        M_CHECK_ARRAY_NON_NULL(tmp_hat_diag, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_loess_weights, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_eval_points, n_genes * 1)
        M_CHECK_ARRAY_NON_NULL(tmp_robust_weights, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_combined_weights, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_residuals, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_indices, n_genes)

        use_quantile_f = use_quantile

        call normalization_pipeline_expert(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            log_transformed_expr = log_transformed_expr,&
            reps_per_tissue = reps_per_tissue,&
            n_tissues = n_tissues,&
            tmp_expr_copy = tmp_expr_copy,&
            tmp_loess_y = tmp_loess_y,&
            tmp_indices_used = tmp_indices_used,&
            tmp_yhat_global = tmp_yhat_global,&
            tmp_int_workspace = tmp_int_workspace,&
            int_workspace_size = int_workspace_size,&
            tmp_real_workspace = tmp_real_workspace,&
            real_workspace_size = real_workspace_size,&
            tmp_hat_diag = tmp_hat_diag,&
            tmp_loess_weights = tmp_loess_weights,&
            tmp_eval_points = tmp_eval_points,&
            tmp_robust_weights = tmp_robust_weights,&
            tmp_combined_weights = tmp_combined_weights,&
            tmp_residuals = tmp_residuals,&
            tmp_permutation_indices = tmp_permutation_indices,&
            span = span,&
            degree = degree,&
            use_quantile = use_quantile_f,&
            ierr = ierr&
        )
    end subroutine normalization_pipeline_expert_c

    !> summary: C-wrapper for [[tox_normalization(module):normalize_by_std_dev(subroutine)]]
    !| This procedure applies a global stabilization based on the relationship between
    !| gene-wise mean expression and empirical standard deviation.
    subroutine normalize_by_std_dev_c(&
            n_genes,&
            n_replicates,&
            expr,&
            normalized_expr,&
            span,&
            degree,&
            ierr&
        ) bind(C, name="normalize_by_std_dev_c")
        use tox_normalization, only: normalize_by_std_dev

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes (rows)
        integer(c_int), intent(in), target :: n_replicates
            !! Number of replicates per gene
        real(c_double), dimension(n_replicates, n_genes), intent(in), target :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_replicates, n_genes), intent(out), target :: normalized_expr
            !! Normalized `expr`
        real(c_double), intent(in), target :: span
            !! LOESS span parameter.
            !! The default value is `0.7_real64`.
        integer(c_int), intent(in), target :: degree
            !! LOESS degree parameter.
            !! The default value is `2_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_replicates)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_ARRAY_NON_NULL(expr, n_replicates * n_genes)
        M_CHECK_ARRAY_NON_NULL(normalized_expr, n_replicates * n_genes)

        call normalize_by_std_dev(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr,&
            span = span,&
            degree = degree,&
            ierr = ierr&
        )
    end subroutine normalize_by_std_dev_c

    !> summary: C-wrapper for [[tox_normalization(module):normalize_by_std_dev_expert(subroutine)]]
    !| This procedure applies a global stabilization based on the relationship between
    !| gene-wise mean expression and empirical standard deviation.
    subroutine normalize_by_std_dev_expert_c(&
            n_genes,&
            n_replicates,&
            expr,&
            normalized_expr,&
            tmp_loess_x,&
            tmp_loess_y,&
            tmp_indices_used,&
            tmp_yhat_global,&
            tmp_int_workspace,&
            int_workspace_size,&
            tmp_real_workspace,&
            real_workspace_size,&
            tmp_hat_diag,&
            tmp_loess_weights,&
            tmp_eval_points,&
            tmp_robust_weights,&
            tmp_combined_weights,&
            tmp_residuals,&
            tmp_permutation_indices,&
            span,&
            degree,&
            ierr&
        ) bind(C, name="normalize_by_std_dev_expert_c")
        use tox_normalization, only: normalize_by_std_dev_expert

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes (rows)
        integer(c_int), intent(in), target :: n_replicates
            !! Number of replicates per gene
        integer(c_int), intent(in), target :: int_workspace_size
            !! Length of integer workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `int_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_genes     |
            !! | save_factorization    | .false.     |
        integer(c_int), intent(in), target :: real_workspace_size
            !! Length of real workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `real_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_genes     |
            !! | save_factorization    | .false.     |
        real(c_double), dimension(n_replicates, n_genes), intent(in), target :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_replicates, n_genes), intent(out), target :: normalized_expr
            !! Normalized `expr`
        real(c_double), dimension(n_genes), intent(out), target :: tmp_loess_x
            !! Work vector for the mean values (X-axis for LOESS)
        real(c_double), dimension(n_genes), intent(out), target :: tmp_loess_y
            !! Work vector for the empirical standard deviations (Y-axis for LOESS)
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_indices_used
            !! Work vector mapping the fitted points back to gene indices
        real(c_double), dimension(n_genes), intent(out), target :: tmp_yhat_global
            !! Work vector for the fitted standard deviations (LOESS predictions)
        integer(c_int), dimension(int_workspace_size), intent(out), target :: tmp_int_workspace
            !! Integer workspace array
        real(c_double), dimension(real_workspace_size), intent(out), target :: tmp_real_workspace
            !! Real workspace array
        real(c_double), dimension(n_genes), intent(out), target :: tmp_hat_diag
            !! Diagonal elements of the LOESS hat matrix
        real(c_double), dimension(n_genes), intent(out), target :: tmp_loess_weights
            !! Per-point weights handed to the LOESS fit
        real(c_double), dimension(n_genes, 1), intent(out), target :: tmp_eval_points
            !! Points the fitted curve is evaluated at
        real(c_double), dimension(n_genes), intent(out), target :: tmp_robust_weights
            !! Robust bisquare weights of the LOESS fit
        real(c_double), dimension(n_genes), intent(out), target :: tmp_combined_weights
            !! Combined weights of the LOESS fit
        real(c_double), dimension(n_genes), intent(out), target :: tmp_residuals
            !! Residuals of the LOESS fit
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_permutation_indices
            !! Permutation indices of the LOESS fit
        real(c_double), intent(in), target :: span
            !! LOESS span parameter.
            !! The default value is `0.7_real64`.
        integer(c_int), intent(in), target :: degree
            !! LOESS degree parameter.
            !! The default value is `2_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_replicates)
        M_CHECK_NON_NULL(int_workspace_size)
        M_CHECK_NON_NULL(real_workspace_size)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_ARRAY_NON_NULL(expr, n_replicates * n_genes)
        M_CHECK_ARRAY_NON_NULL(normalized_expr, n_replicates * n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_loess_x, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_loess_y, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_indices_used, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_yhat_global, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_int_workspace, int_workspace_size)
        M_CHECK_ARRAY_NON_NULL(tmp_real_workspace, real_workspace_size)
        M_CHECK_ARRAY_NON_NULL(tmp_hat_diag, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_loess_weights, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_eval_points, n_genes * 1)
        M_CHECK_ARRAY_NON_NULL(tmp_robust_weights, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_combined_weights, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_residuals, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_indices, n_genes)

        call normalize_by_std_dev_expert(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr,&
            tmp_loess_x = tmp_loess_x,&
            tmp_loess_y = tmp_loess_y,&
            tmp_indices_used = tmp_indices_used,&
            tmp_yhat_global = tmp_yhat_global,&
            tmp_int_workspace = tmp_int_workspace,&
            int_workspace_size = int_workspace_size,&
            tmp_real_workspace = tmp_real_workspace,&
            real_workspace_size = real_workspace_size,&
            tmp_hat_diag = tmp_hat_diag,&
            tmp_loess_weights = tmp_loess_weights,&
            tmp_eval_points = tmp_eval_points,&
            tmp_robust_weights = tmp_robust_weights,&
            tmp_combined_weights = tmp_combined_weights,&
            tmp_residuals = tmp_residuals,&
            tmp_permutation_indices = tmp_permutation_indices,&
            span = span,&
            degree = degree,&
            ierr = ierr&
        )
    end subroutine normalize_by_std_dev_expert_c

    !> summary: C-wrapper for [[tox_normalization(module):root_mean_sq_normalization(subroutine)]]
    !| across tissues (not classical standard deviation).
    subroutine root_mean_sq_normalization_c(&
            n_genes,&
            n_replicates,&
            expr,&
            normalized_expr,&
            ierr&
        ) bind(C, name="root_mean_sq_normalization_c")
        use tox_normalization, only: root_mean_sq_normalization

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes (rows)
        integer(c_int), intent(in), target :: n_replicates
            !! Number of replicates per gene
        real(c_double), dimension(n_replicates, n_genes), intent(in), target :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_replicates, n_genes), intent(out), target :: normalized_expr
            !! Normalized `expr`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_replicates)
        M_CHECK_ARRAY_NON_NULL(expr, n_replicates * n_genes)
        M_CHECK_ARRAY_NON_NULL(normalized_expr, n_replicates * n_genes)

        call root_mean_sq_normalization(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr,&
            ierr = ierr&
        )
    end subroutine root_mean_sq_normalization_c

    !> summary: C-wrapper for [[tox_normalization(module):quantile_normalization(subroutine)]]
    !| Computes average expression per rank across tissues.
    subroutine quantile_normalization_c(&
            n_genes,&
            n_replicates,&
            expr,&
            normalized_expr,&
            rank_means,&
            ierr&
        ) bind(C, name="quantile_normalization_c")
        use tox_normalization, only: quantile_normalization

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes (rows)
        integer(c_int), intent(in), target :: n_replicates
            !! Number of replicates per gene
        real(c_double), dimension(n_replicates, n_genes), intent(in), target :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_replicates, n_genes), intent(out), target :: normalized_expr
            !! Normalized `expr`
        real(c_double), dimension(n_genes), intent(out), target :: rank_means
            !! The mean of each rank across tissues, one per gene
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_replicates)
        M_CHECK_ARRAY_NON_NULL(expr, n_replicates * n_genes)
        M_CHECK_ARRAY_NON_NULL(normalized_expr, n_replicates * n_genes)
        M_CHECK_ARRAY_NON_NULL(rank_means, n_genes)

        call quantile_normalization(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr,&
            rank_means = rank_means,&
            ierr = ierr&
        )
    end subroutine quantile_normalization_c

    !> summary: C-wrapper for [[tox_normalization(module):quantile_normalization_expert(subroutine)]]
    !| Computes average expression per rank across tissues.
    subroutine quantile_normalization_expert_c(&
            n_genes,&
            n_replicates,&
            expr,&
            normalized_expr,&
            rank_means,&
            tmp_genes_row,&
            tmp_perm,&
            ierr&
        ) bind(C, name="quantile_normalization_expert_c")
        use tox_normalization, only: quantile_normalization_expert

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes (rows)
        integer(c_int), intent(in), target :: n_replicates
            !! Number of replicates per gene
        real(c_double), dimension(n_replicates, n_genes), intent(in), target :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_replicates, n_genes), intent(out), target :: normalized_expr
            !! Normalized `expr`
        real(c_double), dimension(n_genes), intent(out), target :: rank_means
            !! The mean of each rank across tissues, one per gene
        real(c_double), dimension(n_genes), intent(out), target :: tmp_genes_row
            !! Temporary vector for sorting a tissue in `expr` across genes
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_perm
            !! Permutation vector
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_replicates)
        M_CHECK_ARRAY_NON_NULL(expr, n_replicates * n_genes)
        M_CHECK_ARRAY_NON_NULL(normalized_expr, n_replicates * n_genes)
        M_CHECK_ARRAY_NON_NULL(rank_means, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_genes_row, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_perm, n_genes)

        call quantile_normalization_expert(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr,&
            rank_means = rank_means,&
            tmp_genes_row = tmp_genes_row,&
            tmp_perm = tmp_perm,&
            ierr = ierr&
        )
    end subroutine quantile_normalization_expert_c

    !> summary: C-wrapper for [[tox_normalization(module):log2_transformation(subroutine)]]
    !| This subroutine performs element-wise `log2(x + 1)` transformation on a
    !| matrix flattened in column-major order. The `log2` is computed via:
    !| `log(x + 1) / log(2)`, which is numerically equivalent and avoids the
    !| non-portable `log2` intrinsic for compatibility with WebAssembly (WASM).
    subroutine log2_transformation_c(&
            n_genes,&
            n_tissues,&
            expr,&
            transformed_expr,&
            ierr&
        ) bind(C, name="log2_transformation_c")
        use tox_normalization, only: log2_transformation

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes (rows)
        integer(c_int), intent(in), target :: n_tissues
            !! Number of tissues
        real(c_double), dimension(n_tissues, n_genes), intent(in), target :: expr
            !! Gene Expression matrix, from [[tox_normalization(module):calc_tiss_avg(subroutine)]]
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_tissues, n_genes), intent(out), target :: transformed_expr
            !! Log-transformed `expr`
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_tissues)
        M_CHECK_ARRAY_NON_NULL(expr, n_tissues * n_genes)
        M_CHECK_ARRAY_NON_NULL(transformed_expr, n_tissues * n_genes)

        call log2_transformation(&
            n_genes = n_genes,&
            n_tissues = n_tissues,&
            expr = expr,&
            transformed_expr = transformed_expr,&
            ierr = ierr&
        )
    end subroutine log2_transformation_c

    !> summary: C-wrapper for [[tox_normalization(module):calc_tiss_avg(subroutine)]]
    !| For each tissue of tissue replicates, this subroutine computes the average
    !| expression per gene.
    subroutine calc_tiss_avg_c(&
            n_genes,&
            n_tissues,&
            reps_per_tissue,&
            expr,&
            tissue_averages,&
            ierr&
        ) bind(C, name="calc_tiss_avg_c")
        use tox_normalization, only: calc_tiss_avg

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes (rows)
        integer(c_int), intent(in), target :: n_tissues
            !! Number of tissues
        integer(c_int), dimension(n_tissues), intent(in), target :: reps_per_tissue
            !! Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
            !! e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(sum(reps_per_tissue), n_genes), intent(in), target :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_tissues, n_genes), intent(out), target :: tissue_averages
            !! Tissue averages per gene
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_tissues)
        M_CHECK_ARRAY_NON_NULL(reps_per_tissue, n_tissues)
        M_CHECK_ARRAY_NON_NULL(expr, (sum(reps_per_tissue)) * n_genes)
        M_CHECK_ARRAY_NON_NULL(tissue_averages, n_tissues * n_genes)

        call calc_tiss_avg(&
            n_genes = n_genes,&
            n_tissues = n_tissues,&
            reps_per_tissue = reps_per_tissue,&
            expr = expr,&
            tissue_averages = tissue_averages,&
            ierr = ierr&
        )
    end subroutine calc_tiss_avg_c

    !> summary: C-wrapper for [[tox_normalization(module):calc_fchange(subroutine)]]
    !| For each control-condition pair, this subroutine computes the `log2 fold change`
    !| by subtracting the expression value in the control group from the corresponding
    !| value in the condition group, for all genes.
    subroutine calc_fchange_c(&
            n_genes,&
            n_tissues,&
            n_pairs,&
            control_tissues,&
            condition_tissues,&
            expr,&
            fold_changes,&
            ierr&
        ) bind(C, name="calc_fchange_c")
        use tox_normalization, only: calc_fchange

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes (rows)
        integer(c_int), intent(in), target :: n_tissues
            !! Number of tissues
        integer(c_int), intent(in), target :: n_pairs
            !! Number of control-condition pairs
        integer(c_int), dimension(n_pairs), intent(in), target :: control_tissues
            !! Control tissue indices
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_tissues`.
        integer(c_int), dimension(n_pairs), intent(in), target :: condition_tissues
            !! Condition tissue indices
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_tissues`.
        real(c_double), dimension(n_tissues, n_genes), intent(in), target :: expr
            !! Gene Expression matrix, from [[tox_normalization(module):calc_tiss_avg(subroutine)]]
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_pairs, n_genes), intent(out), target :: fold_changes
            !! Output matrix for fold changes
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_tissues)
        M_CHECK_NON_NULL(n_pairs)
        M_CHECK_ARRAY_NON_NULL(control_tissues, n_pairs)
        M_CHECK_ARRAY_NON_NULL(condition_tissues, n_pairs)
        M_CHECK_ARRAY_NON_NULL(expr, n_tissues * n_genes)
        M_CHECK_ARRAY_NON_NULL(fold_changes, n_pairs * n_genes)

        call calc_fchange(&
            n_genes = n_genes,&
            n_tissues = n_tissues,&
            n_pairs = n_pairs,&
            control_tissues = control_tissues,&
            condition_tissues = condition_tissues,&
            expr = expr,&
            fold_changes = fold_changes,&
            ierr = ierr&
        )
    end subroutine calc_fchange_c

end module tox_normalization_c
#endif
