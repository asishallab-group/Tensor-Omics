#include <src/macros.h>

!> summary: Wrappers for [[tox_normalization_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_normalization
    use tox_normalization_kernel, only: calc_fchange_kernel, calc_tiss_avg_kernel, log2_transformation_kernel, normalization_pipeline_alloc_kernel
    use tox_normalization_kernel, only: normalize_by_std_dev_alloc_kernel, normalize_unit_length_kernel, quantile_normalization_kernel, root_mean_sq_normalization_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, clear_err_arg_pos
    use tox_errors, only: set_err, validate_all_in_range_int, validate_dimension_size, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: normalize_unit_length
    public :: normalization_pipeline_alloc
    public :: normalize_by_std_dev_alloc
    public :: root_mean_sq_normalization
    public :: quantile_normalization
    public :: quantile_normalization_alloc
    public :: log2_transformation
    public :: calc_tiss_avg
    public :: calc_fchange

contains

    !> summary: Validates its inputs, then calls [[tox_normalization_kernel(module):normalize_unit_length_kernel]].
    subroutine normalize_unit_length(&
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
        call validate_dimension_size(n_dims, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return

        call normalize_unit_length_kernel(&
            vector = vector,&
            n_dims = n_dims,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine normalize_unit_length

    !> summary: Validates its inputs, then calls [[tox_normalization_kernel(module):normalization_pipeline_alloc_kernel]].
    !| Final result is in log_transformed_expr. If fold change is needed, call calc_fchange separately.
    subroutine normalization_pipeline_alloc(&
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
            !! e.g. `[2,3]` means `5` total replicates per gene, with the `expr(1:2, i_gene)` related to the first tissue and `expr(3:, i_gene)` related to the second one.
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
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_replicates, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_tissues, ierr, arg_pos=6_int32)
        call validate_in_range_real(span, ierr, arg_pos=7_int32)
        if (is_err(ierr)) return

        call normalization_pipeline_alloc_kernel(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            log_transformed_expr = log_transformed_expr,&
            reps_per_tissue = reps_per_tissue,&
            n_tissues = n_tissues,&
            span = span,&
            degree = degree,&
            use_quantile = use_quantile,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine normalization_pipeline_alloc

    !> summary: Validates its inputs, then calls [[tox_normalization_kernel(module):normalize_by_std_dev_alloc_kernel]].
    !| This procedure applies a global stabilization based on the relationship between
    !| gene-wise mean expression and empirical standard deviation.
    subroutine normalize_by_std_dev_alloc(&
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

        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_replicates, ierr, arg_pos=2_int32)
        call validate_in_range_real(span, ierr, arg_pos=5_int32)
        if (is_err(ierr)) return

        call normalize_by_std_dev_alloc_kernel(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr,&
            span = span,&
            degree = degree,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine normalize_by_std_dev_alloc

    !> summary: Validates its inputs, then calls [[tox_normalization_kernel(module):root_mean_sq_normalization_kernel]].
    !| across tissues (not classical standard deviation).
    subroutine root_mean_sq_normalization(&
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
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_replicates, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return

        call root_mean_sq_normalization_kernel(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr&
        )
    end subroutine root_mean_sq_normalization

    !> summary: Validates its inputs, then calls [[tox_normalization_kernel(module):quantile_normalization_kernel]].
    !| Computes average expression per rank across tissues.
    subroutine quantile_normalization(&
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
            !! Preallocated vector to store rank means
        real(real64), dimension(n_genes), intent(out) :: tmp_genes_row
            !! Temporary vector for sorting a tissue in `expr` across genes
        integer(int32), dimension(n_genes), intent(out) :: tmp_perm
            !! Permutation vector
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_replicates, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return

        call quantile_normalization_kernel(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr,&
            rank_means = rank_means,&
            tmp_genes_row = tmp_genes_row,&
            tmp_perm = tmp_perm&
        )
    end subroutine quantile_normalization

    !> summary: Allocates its work arrays, then calls [[tox_normalization_kernel(module):quantile_normalization_kernel]].
    !| Computes average expression per rank across tissues.
    subroutine quantile_normalization_alloc(&
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
            !! Preallocated vector to store rank means
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        real(real64), dimension(:), allocatable :: tmp_genes_row
        integer(int32), dimension(:), allocatable :: tmp_perm

        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_replicates, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_genes_row(n_genes))
        M_ALLOCATE(tmp_perm(n_genes))

        call quantile_normalization_kernel(&
            n_genes = n_genes,&
            n_replicates = n_replicates,&
            expr = expr,&
            normalized_expr = normalized_expr,&
            rank_means = rank_means,&
            tmp_genes_row = tmp_genes_row,&
            tmp_perm = tmp_perm&
        )
    end subroutine quantile_normalization_alloc

    !> summary: Validates its inputs, then calls [[tox_normalization_kernel(module):log2_transformation_kernel]].
    !| This subroutine performs element-wise `log2(x + 1)` transformation on a
    !| matrix flattened in column-major order. The `log2` is computed via:
    !| `log(x + 1) / log(2)`, which is numerically equivalent and avoids the
    !| non-portable `log2` intrinsic for compatibility with WebAssembly (WASM).
    subroutine log2_transformation(&
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
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_tissues, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return

        call log2_transformation_kernel(&
            n_genes = n_genes,&
            n_tissues = n_tissues,&
            expr = expr,&
            transformed_expr = transformed_expr,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine log2_transformation

    !> summary: Validates its inputs, then calls [[tox_normalization_kernel(module):calc_tiss_avg_kernel]].
    !| For each tissue of tissue replicates, this subroutine computes the average
    !| expression per gene.
    subroutine calc_tiss_avg(&
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
            !! e.g. `[2,3]` means `5` total replicates per gene, with the `expr(1:2, i_gene)` related to the first tissue and `expr(3:, i_gene)` related to the second one.
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
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_tissues, ierr, arg_pos=2_int32)
        call validate_all_in_range_int(reps_per_tissue, n_tissues, ierr, arg_pos=3_int32, min=1_int32)
        if (is_err(ierr)) return

        call calc_tiss_avg_kernel(&
            n_genes = n_genes,&
            n_tissues = n_tissues,&
            reps_per_tissue = reps_per_tissue,&
            expr = expr,&
            tissue_averages = tissue_averages&
        )
    end subroutine calc_tiss_avg

    !> summary: Validates its inputs, then calls [[tox_normalization_kernel(module):calc_fchange_kernel]].
    !| For each control-condition pair, this subroutine computes the `log2 fold change`
    !| by subtracting the expression value in the control group from the corresponding
    !| value in the condition group, for all genes.
    subroutine calc_fchange(&
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
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_tissues, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_pairs, ierr, arg_pos=3_int32)
        call validate_all_in_range_int(control_tissues, n_pairs, ierr, arg_pos=4_int32, min=1_int32, max=n_tissues)
        call validate_all_in_range_int(condition_tissues, n_pairs, ierr, arg_pos=5_int32, min=1_int32, max=n_tissues)
        if (is_err(ierr)) return

        call calc_fchange_kernel(&
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
