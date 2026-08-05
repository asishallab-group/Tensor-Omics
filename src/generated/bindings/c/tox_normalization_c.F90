#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_normalization(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_normalization_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: normalize_unit_length_c
    public :: normalization_pipeline_c
    public :: normalize_by_std_dev_c
    public :: root_mean_sq_normalization_c
    public :: quantile_normalization_expert_c
    public :: quantile_normalization_c
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

    !> summary: C-wrapper for [[tox_normalization(module):normalization_pipeline_alloc(subroutine)]]
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
        use tox_normalization, only: normalization_pipeline_alloc

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
            !! e.g. `[2,3]` means `5` total replicates per gene, with the `expr(1:2, i_gene)` related to the first tissue and `expr(3:, i_gene)` related to the second one.
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

        call normalization_pipeline_alloc(&
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

    !> summary: C-wrapper for [[tox_normalization(module):normalize_by_std_dev_alloc(subroutine)]]
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
        use tox_normalization, only: normalize_by_std_dev_alloc

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

        call normalize_by_std_dev_alloc(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr,&
            span = span,&
            degree = degree,&
            ierr = ierr&
        )
    end subroutine normalize_by_std_dev_c

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
            !! Preallocated vector to store rank means
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

        call quantile_normalization(&
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

    !> summary: C-wrapper for [[tox_normalization(module):quantile_normalization_alloc(subroutine)]]
    !| Computes average expression per rank across tissues.
    subroutine quantile_normalization_c(&
            n_genes,&
            n_replicates,&
            expr,&
            normalized_expr,&
            rank_means,&
            ierr&
        ) bind(C, name="quantile_normalization_c")
        use tox_normalization, only: quantile_normalization_alloc

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
            !! Preallocated vector to store rank means
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_replicates)
        M_CHECK_ARRAY_NON_NULL(expr, n_replicates * n_genes)
        M_CHECK_ARRAY_NON_NULL(normalized_expr, n_replicates * n_genes)
        M_CHECK_ARRAY_NON_NULL(rank_means, n_genes)

        call quantile_normalization_alloc(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr,&
            rank_means = rank_means,&
            ierr = ierr&
        )
    end subroutine quantile_normalization_c

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
            !! e.g. `[2,3]` means `5` total replicates per gene, with the `expr(1:2, i_gene)` related to the first tissue and `expr(3:, i_gene)` related to the second one.
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
