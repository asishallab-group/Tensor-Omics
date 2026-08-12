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
!|
!| Generated from [[tox_normalization_impl(module)]]; do not edit -- regenerate instead.
module tox_normalization
    use tox_normalization_impl, only: calc_fchange_impl, calc_tiss_avg_impl, log2_transformation_impl, normalization_pipeline_impl
    use tox_normalization_impl, only: normalize_by_std_dev_impl, normalize_unit_length_impl, quantile_normalization_impl, root_mean_sq_normalization_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_loess_impl, only: tox_loess_required_workspace
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, clear_err_arg_pos
    use tox_errors, only: set_err, validate_all_in_range_int, validate_dimension_size, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: normalize_unit_length
    public :: normalization_pipeline
    public :: normalization_pipeline_expert
    public :: normalize_by_std_dev
    public :: normalize_by_std_dev_expert
    public :: root_mean_sq_normalization
    public :: quantile_normalization
    public :: quantile_normalization_expert
    public :: log2_transformation
    public :: calc_tiss_avg
    public :: calc_fchange

contains

    !> summary: Validates its inputs, then calls [[tox_normalization_impl(module):normalize_unit_length_impl]].
    pure subroutine normalize_unit_length(&
            vector,&
            n_dims,&
            ierr&
        )
        integer(int32), intent(in) :: n_dims
            !! number of elements in `vector`
        real(real64), dimension(n_dims), intent(inout) :: vector
            !! Vector that will be normalized to unit length
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dims, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return
#endif

        call normalize_unit_length_impl(&
            vector = vector,&
            n_dims = n_dims,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine normalize_unit_length

    !> summary: Validates its inputs, prepares what [[tox_normalization_impl(module):normalization_pipeline_impl]] needs, then calls it. The entry point to reach for first; see [[tox_normalization(module):normalization_pipeline_expert]] to prepare it yourself.
    !| Final result is in log_transformed_expr. If fold change is needed, call calc_fchange separately.
    subroutine normalization_pipeline(&
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
        )
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_tissues, n_genes), intent(out) :: log_transformed_expr
            !! Log-transformed grouped `expr`
        integer(int32), dimension(n_tissues), intent(in) :: reps_per_tissue
            !! Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
            !! e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
        real(real64), intent(in), optional :: span
            !! LOESS span parameter.
            !! The default value is `0.7_real64`.
        integer(int32), intent(in), optional :: degree
            !! LOESS degree parameter.
            !! The default value is `2_int32`.
        logical, intent(in), optional :: use_quantile
            !! Use quantile normalization.
            !! The default value is `.false.`.
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), dimension(:, :), allocatable :: tmp_expr_copy
        real(real64), dimension(:), allocatable :: tmp_loess_y
        integer(int32), dimension(:), allocatable :: tmp_indices_used
        real(real64), dimension(:), allocatable :: tmp_yhat_global
        integer(int32), dimension(:), allocatable :: tmp_int_workspace
        integer(int32) :: int_workspace_size
        real(real64), dimension(:), allocatable :: tmp_real_workspace
        integer(int32) :: real_workspace_size
        real(real64), dimension(:), allocatable :: tmp_hat_diag
        real(real64), dimension(:), allocatable :: tmp_loess_weights
        real(real64), dimension(:, :), allocatable :: tmp_eval_points
        real(real64), dimension(:), allocatable :: tmp_robust_weights
        real(real64), dimension(:), allocatable :: tmp_combined_weights
        real(real64), dimension(:), allocatable :: tmp_residuals
        integer(int32), dimension(:), allocatable :: tmp_permutation_indices

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_replicates, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_tissues, ierr, arg_pos=6_int32)
        call validate_in_range_real(span, ierr, arg_pos=7_int32)
        if (is_err(ierr)) return
#endif

        call tox_loess_required_workspace(&
            n_dim = 1_int32,&
            max_neighborhood_size = n_genes,&
            int_workspace_size = int_workspace_size,&
            real_workspace_size = real_workspace_size,&
            save_factorization = .false.&
        )
        M_ALLOCATE(tmp_expr_copy(n_replicates, n_genes))
        M_ALLOCATE(tmp_loess_y(n_genes))
        M_ALLOCATE(tmp_indices_used(n_genes))
        M_ALLOCATE(tmp_yhat_global(n_genes))
        M_ALLOCATE(tmp_int_workspace(int_workspace_size))
        M_ALLOCATE(tmp_real_workspace(real_workspace_size))
        M_ALLOCATE(tmp_hat_diag(n_genes))
        M_ALLOCATE(tmp_loess_weights(n_genes))
        M_ALLOCATE(tmp_eval_points(n_genes, 1))
        M_ALLOCATE(tmp_robust_weights(n_genes))
        M_ALLOCATE(tmp_combined_weights(n_genes))
        M_ALLOCATE(tmp_residuals(n_genes))
        M_ALLOCATE(tmp_permutation_indices(n_genes))

        call normalization_pipeline_impl(&
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
            use_quantile = use_quantile,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine normalization_pipeline

    !> summary: Validates its inputs, then calls [[tox_normalization_impl(module):normalization_pipeline_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_normalization(module):normalization_pipeline]] does both.
    !| Final result is in log_transformed_expr. If fold change is needed, call calc_fchange separately.
    subroutine normalization_pipeline_expert(&
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
        )
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        integer(int32), intent(in) :: int_workspace_size
            !! Length of integer workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `int_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_genes     |
            !! | save_factorization    | .false.     |
        integer(int32), intent(in) :: real_workspace_size
            !! Length of real workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `real_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_genes     |
            !! | save_factorization    | .false.     |
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_tissues, n_genes), intent(out) :: log_transformed_expr
            !! Log-transformed grouped `expr`
        integer(int32), dimension(n_tissues), intent(in) :: reps_per_tissue
            !! Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
            !! e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
        real(real64), dimension(n_replicates, n_genes), intent(out) :: tmp_expr_copy
            !! Work matrix the pipeline normalizes in place
        real(real64), dimension(n_genes), intent(out) :: tmp_loess_y
            !! Work vector for the empirical standard deviations (Y-axis for LOESS)
        integer(int32), dimension(n_genes), intent(out) :: tmp_indices_used
            !! Work vector mapping the fitted points back to gene indices
        real(real64), dimension(n_genes), intent(out) :: tmp_yhat_global
            !! Work vector for the fitted standard deviations (LOESS predictions)
        integer(int32), dimension(int_workspace_size), intent(out) :: tmp_int_workspace
            !! Integer workspace array
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
            !! The default value is `0.7_real64`.
        integer(int32), intent(in), optional :: degree
            !! LOESS degree parameter.
            !! The default value is `2_int32`.
        logical, intent(in), optional :: use_quantile
            !! Use quantile normalization.
            !! The default value is `.false.`.
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_replicates, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_tissues, ierr, arg_pos=6_int32)
        call validate_dimension_size(int_workspace_size, ierr, arg_pos=12_int32)
        call validate_dimension_size(real_workspace_size, ierr, arg_pos=14_int32)
        call validate_in_range_real(span, ierr, arg_pos=22_int32)
        if (is_err(ierr)) return
#endif

        call normalization_pipeline_impl(&
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
            use_quantile = use_quantile,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine normalization_pipeline_expert

    !> summary: Validates its inputs, prepares what [[tox_normalization_impl(module):normalize_by_std_dev_impl]] needs, then calls it. The entry point to reach for first; see [[tox_normalization(module):normalize_by_std_dev_expert]] to prepare it yourself.
    !| This procedure applies a global stabilization based on the relationship between
    !| gene-wise mean expression and empirical standard deviation.
    subroutine normalize_by_std_dev(&
            n_genes,&
            n_replicates,&
            expr,&
            normalized_expr,&
            span,&
            degree,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_replicates, n_genes), intent(out) :: normalized_expr
            !! Normalized `expr`
        real(real64), intent(in), optional :: span
            !! LOESS span parameter.
            !! The default value is `0.7_real64`.
        integer(int32), intent(in), optional :: degree
            !! LOESS degree parameter.
            !! The default value is `2_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), dimension(:), allocatable :: tmp_loess_x
        real(real64), dimension(:), allocatable :: tmp_loess_y
        integer(int32), dimension(:), allocatable :: tmp_indices_used
        real(real64), dimension(:), allocatable :: tmp_yhat_global
        integer(int32), dimension(:), allocatable :: tmp_int_workspace
        integer(int32) :: int_workspace_size
        real(real64), dimension(:), allocatable :: tmp_real_workspace
        integer(int32) :: real_workspace_size
        real(real64), dimension(:), allocatable :: tmp_hat_diag
        real(real64), dimension(:), allocatable :: tmp_loess_weights
        real(real64), dimension(:, :), allocatable :: tmp_eval_points
        real(real64), dimension(:), allocatable :: tmp_robust_weights
        real(real64), dimension(:), allocatable :: tmp_combined_weights
        real(real64), dimension(:), allocatable :: tmp_residuals
        integer(int32), dimension(:), allocatable :: tmp_permutation_indices

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_replicates, ierr, arg_pos=2_int32)
        call validate_in_range_real(span, ierr, arg_pos=5_int32)
        if (is_err(ierr)) return
#endif

        call tox_loess_required_workspace(&
            n_dim = 1_int32,&
            max_neighborhood_size = n_genes,&
            int_workspace_size = int_workspace_size,&
            real_workspace_size = real_workspace_size,&
            save_factorization = .false.&
        )
        M_ALLOCATE(tmp_loess_x(n_genes))
        M_ALLOCATE(tmp_loess_y(n_genes))
        M_ALLOCATE(tmp_indices_used(n_genes))
        M_ALLOCATE(tmp_yhat_global(n_genes))
        M_ALLOCATE(tmp_int_workspace(int_workspace_size))
        M_ALLOCATE(tmp_real_workspace(real_workspace_size))
        M_ALLOCATE(tmp_hat_diag(n_genes))
        M_ALLOCATE(tmp_loess_weights(n_genes))
        M_ALLOCATE(tmp_eval_points(n_genes, 1))
        M_ALLOCATE(tmp_robust_weights(n_genes))
        M_ALLOCATE(tmp_combined_weights(n_genes))
        M_ALLOCATE(tmp_residuals(n_genes))
        M_ALLOCATE(tmp_permutation_indices(n_genes))

        call normalize_by_std_dev_impl(&
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
        call clear_err_arg_pos(ierr)
    end subroutine normalize_by_std_dev

    !> summary: Validates its inputs, then calls [[tox_normalization_impl(module):normalize_by_std_dev_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_normalization(module):normalize_by_std_dev]] does both.
    !| This procedure applies a global stabilization based on the relationship between
    !| gene-wise mean expression and empirical standard deviation.
    subroutine normalize_by_std_dev_expert(&
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
        )
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        integer(int32), intent(in) :: int_workspace_size
            !! Length of integer workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `int_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_genes     |
            !! | save_factorization    | .false.     |
        integer(int32), intent(in) :: real_workspace_size
            !! Length of real workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `real_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_genes     |
            !! | save_factorization    | .false.     |
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_replicates, n_genes), intent(out) :: normalized_expr
            !! Normalized `expr`
        real(real64), dimension(n_genes), intent(out) :: tmp_loess_x
            !! Work vector for the mean values (X-axis for LOESS)
        real(real64), dimension(n_genes), intent(out) :: tmp_loess_y
            !! Work vector for the empirical standard deviations (Y-axis for LOESS)
        integer(int32), dimension(n_genes), intent(out) :: tmp_indices_used
            !! Work vector mapping the fitted points back to gene indices
        real(real64), dimension(n_genes), intent(out) :: tmp_yhat_global
            !! Work vector for the fitted standard deviations (LOESS predictions)
        integer(int32), dimension(int_workspace_size), intent(out) :: tmp_int_workspace
            !! Integer workspace array
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
            !! The default value is `0.7_real64`.
        integer(int32), intent(in), optional :: degree
            !! LOESS degree parameter.
            !! The default value is `2_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_replicates, ierr, arg_pos=2_int32)
        call validate_dimension_size(int_workspace_size, ierr, arg_pos=10_int32)
        call validate_dimension_size(real_workspace_size, ierr, arg_pos=12_int32)
        call validate_in_range_real(span, ierr, arg_pos=20_int32)
        if (is_err(ierr)) return
#endif

        call normalize_by_std_dev_impl(&
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
        call clear_err_arg_pos(ierr)
    end subroutine normalize_by_std_dev_expert

    !> summary: Validates its inputs, then calls [[tox_normalization_impl(module):root_mean_sq_normalization_impl]].
    !| across tissues (not classical standard deviation).
    pure subroutine root_mean_sq_normalization(&
            n_genes,&
            n_replicates,&
            expr,&
            normalized_expr,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_replicates, n_genes), intent(out) :: normalized_expr
            !! Normalized `expr`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_replicates, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return
#endif

        call root_mean_sq_normalization_impl(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr&
        )
    end subroutine root_mean_sq_normalization

    !> summary: Validates its inputs, prepares what [[tox_normalization_impl(module):quantile_normalization_impl]] needs, then calls it. The entry point to reach for first; see [[tox_normalization(module):quantile_normalization_expert]] to prepare it yourself.
    !| Computes average expression per rank across tissues.
    pure subroutine quantile_normalization(&
            n_genes,&
            n_replicates,&
            expr,&
            normalized_expr,&
            rank_means,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_replicates, n_genes), intent(out) :: normalized_expr
            !! Normalized `expr`
        real(real64), dimension(n_genes), intent(out) :: rank_means
            !! The mean of each rank across tissues, one per gene
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        real(real64), dimension(:), allocatable :: tmp_genes_row
        integer(int32), dimension(:), allocatable :: tmp_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_replicates, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_genes_row(n_genes))
        M_ALLOCATE(tmp_perm(n_genes))

        call quantile_normalization_impl(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr,&
            rank_means = rank_means,&
            tmp_genes_row = tmp_genes_row,&
            tmp_perm = tmp_perm&
        )
    end subroutine quantile_normalization

    !> summary: Validates its inputs, then calls [[tox_normalization_impl(module):quantile_normalization_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_normalization(module):quantile_normalization]] does both.
    !| Computes average expression per rank across tissues.
    pure subroutine quantile_normalization_expert(&
            n_genes,&
            n_replicates,&
            expr,&
            normalized_expr,&
            rank_means,&
            tmp_genes_row,&
            tmp_perm,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene
        real(real64), dimension(n_replicates, n_genes), intent(in) :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_replicates, n_genes), intent(out) :: normalized_expr
            !! Normalized `expr`
        real(real64), dimension(n_genes), intent(out) :: rank_means
            !! The mean of each rank across tissues, one per gene
        real(real64), dimension(n_genes), intent(out) :: tmp_genes_row
            !! Temporary vector for sorting a tissue in `expr` across genes
        integer(int32), dimension(n_genes), intent(out) :: tmp_perm
            !! Permutation vector
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_replicates, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return
#endif

        call quantile_normalization_impl(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr,&
            rank_means = rank_means,&
            tmp_genes_row = tmp_genes_row,&
            tmp_perm = tmp_perm&
        )
    end subroutine quantile_normalization_expert

    !> summary: Validates its inputs, then calls [[tox_normalization_impl(module):log2_transformation_impl]].
    !| This subroutine performs element-wise `log2(x + 1)` transformation on a
    !| matrix flattened in column-major order. The `log2` is computed via:
    !| `log(x + 1) / log(2)`, which is numerically equivalent and avoids the
    !| non-portable `log2` intrinsic for compatibility with WebAssembly (WASM).
    pure subroutine log2_transformation(&
            n_genes,&
            n_tissues,&
            expr,&
            transformed_expr,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        real(real64), dimension(n_tissues, n_genes), intent(in) :: expr
            !! Gene Expression matrix, from [[tox_normalization(module):calc_tiss_avg(subroutine)]]
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_tissues, n_genes), intent(out) :: transformed_expr
            !! Log-transformed `expr`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_tissues, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return
#endif

        call log2_transformation_impl(&
            n_genes = n_genes,&
            n_tissues = n_tissues,&
            expr = expr,&
            transformed_expr = transformed_expr,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine log2_transformation

    !> summary: Validates its inputs, then calls [[tox_normalization_impl(module):calc_tiss_avg_impl]].
    !| For each tissue of tissue replicates, this subroutine computes the average
    !| expression per gene.
    pure subroutine calc_tiss_avg(&
            n_genes,&
            n_tissues,&
            reps_per_tissue,&
            expr,&
            tissue_averages,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        integer(int32), dimension(n_tissues), intent(in) :: reps_per_tissue
            !! Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
            !! e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(sum(reps_per_tissue), n_genes), intent(in) :: expr
            !! Gene Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_tissues, n_genes), intent(out) :: tissue_averages
            !! Tissue averages per gene
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_tissues, ierr, arg_pos=2_int32)
        call validate_all_in_range_int(reps_per_tissue, n_tissues, ierr, arg_pos=3_int32, min=1_int32)
        if (is_err(ierr)) return
#endif

        call calc_tiss_avg_impl(&
            n_genes = n_genes,&
            n_tissues = n_tissues,&
            reps_per_tissue = reps_per_tissue,&
            expr = expr,&
            tissue_averages = tissue_averages&
        )
    end subroutine calc_tiss_avg

    !> summary: Validates its inputs, then calls [[tox_normalization_impl(module):calc_fchange_impl]].
    !| For each control-condition pair, this subroutine computes the `log2 fold change`
    !| by subtracting the expression value in the control group from the corresponding
    !| value in the condition group, for all genes.
    pure subroutine calc_fchange(&
            n_genes,&
            n_tissues,&
            n_pairs,&
            control_tissues,&
            condition_tissues,&
            expr,&
            fold_changes,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! Number of genes (rows)
        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        integer(int32), intent(in) :: n_pairs
            !! Number of control-condition pairs
        integer(int32), dimension(n_pairs), intent(in) :: control_tissues
            !! Control tissue indices
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_tissues`.
        integer(int32), dimension(n_pairs), intent(in) :: condition_tissues
            !! Condition tissue indices
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_tissues`.
        real(real64), dimension(n_tissues, n_genes), intent(in) :: expr
            !! Gene Expression matrix, from [[tox_normalization(module):calc_tiss_avg(subroutine)]]
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_pairs, n_genes), intent(out) :: fold_changes
            !! Output matrix for fold changes
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_tissues, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_pairs, ierr, arg_pos=3_int32)
        call validate_all_in_range_int(control_tissues, n_pairs, ierr, arg_pos=4_int32, min=1_int32, max=n_tissues)
        call validate_all_in_range_int(condition_tissues, n_pairs, ierr, arg_pos=5_int32, min=1_int32, max=n_tissues)
        if (is_err(ierr)) return
#endif

        call calc_fchange_impl(&
            n_genes = n_genes,&
            n_tissues = n_tissues,&
            n_pairs = n_pairs,&
            control_tissues = control_tissues,&
            condition_tissues = condition_tissues,&
            expr = expr,&
            fold_changes = fold_changes&
        )
    end subroutine calc_fchange

end module tox_normalization
