#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_get_outliers(module)]]
!| Gene outliers, from how far each gene sits from its family's centroid.
!|
!| The pipeline is three steps, each callable on its own. `compute_rdi` turns raw distances into
!| a relative distance index, scaled per family so families of different spread are comparable.
!| `compute_family_scaling` fits that scaling with LOESS against family size. `identify_outliers`
!| applies the threshold and reports which genes exceed it.
!|
!| `detect_outliers` runs all three in one call, and is the entry point to reach for first.
module tox_get_outliers_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL, ERR_INVALID_INPUT, ERR_ALLOC_FAIL
    M_IMPLICIT_NONE
    private

    public :: compute_family_scaling_c
    public :: compute_family_scaling_expert_c
    public :: compute_rdi_c
    public :: compute_rdi_expert_c
    public :: identify_outliers_c
    public :: detect_outliers_c
    public :: detect_outliers_expert_c

contains

    !> summary: C-wrapper for [[tox_get_outliers(module):compute_family_scaling(subroutine)]]
    !| Uses LOESS on the median/stddev of intra-family distances for scaling, regardless of orthologs.
    subroutine compute_family_scaling_c(&
            n_genes,&
            n_families,&
            distances,&
            gene_to_fam,&
            dscale,&
            loess_x,&
            loess_y,&
            indices_used,&
            span,&
            degree,&
            mode,&
            n_iters,&
            low_sd_cutoff,&
            excluded_low_sd,&
            ierr&
        ) bind(C, name="compute_family_scaling_c")
        use tox_get_outliers, only: compute_family_scaling
        use tox_loess_impl, only: MODE_PLAIN, MODE_ROBUST

        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes
        integer(c_int), intent(in), target :: n_families
            !! Total number of gene families
        real(c_double), dimension(n_genes), intent(in), target :: distances
            !! Array of Euclidean distances for each gene
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Mapping of each gene to its family (1-based)
        real(c_double), dimension(n_families), intent(out), target :: dscale
            !! Array of scaling factors per family (output)
        real(c_double), dimension(n_families), intent(out), target :: loess_x
            !! Reference x-coordinates for LOESS smoothing
        real(c_double), dimension(n_families), intent(out), target :: loess_y
            !! Reference y-coordinates for LOESS smoothing
        integer(c_int), dimension(n_families), intent(out), target :: indices_used
            !! Indices of reference points used for smoothing
        real(c_double), intent(in), target :: span
            !! Span parameter for LOESS smoothing, passed straight to
            !! [[tox_loess_impl(module):loess_fit_plain_impl(subroutine)]], so it is held to that
            !! procedure's own range rather than to the NaN tolerance the distance data carries.
            !! The default value is `0.7_real64`.
            !! The minimum valid value is `EPS_LOESS`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), intent(in), target :: degree
            !! Degree of the LOESS polynomial
            !! The default value is `2_int32`.
        character(len=1, kind=c_char), dimension(6), intent(in), target :: mode
            !! Mode for LOESS fitting
            !! The default value is `'robust'`.
            !!
            !! | Mode                 | Value                                            |
            !! |----------------------|--------------------------------------------------|
            !! | Plain LOESS fitting  | [[tox_loess_impl(module):MODE_PLAIN(variable)]]  |
            !! | Robust LOESS fitting | [[tox_loess_impl(module):MODE_ROBUST(variable)]] |
        integer(c_int), intent(in), target :: n_iters
            !! Number of iterations for robust LOESS fitting
            !! The default value is `3_int32`.
        real(c_double), intent(out), target :: low_sd_cutoff
            !! cutoff used to filter families with low std
        integer(c_int), dimension(n_families), intent(out), target :: excluded_low_sd
            !! Mask to save those families that have low sd
        integer(c_int), intent(out), target :: ierr
            !! Error code
        integer(int32) :: mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_NON_NULL(n_iters)
        M_CHECK_NON_NULL(low_sd_cutoff)
        M_CHECK_ARRAY_NON_NULL(distances, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(dscale, n_families)
        M_CHECK_ARRAY_NON_NULL(loess_x, n_families)
        M_CHECK_ARRAY_NON_NULL(loess_y, n_families)
        M_CHECK_ARRAY_NON_NULL(indices_used, n_families)
        M_CHECK_ARRAY_NON_NULL(mode, 6)
        M_CHECK_ARRAY_NON_NULL(excluded_low_sd, n_families)

        block
            character(len=:), allocatable :: mode_f
            call c_char_1d_as_string(mode, mode_f, ierr)
            if (is_err(ierr)) return

            select case (mode_f)
                case ("plain")
                    mode_mode_f = MODE_PLAIN
                case ("robust")
                    mode_mode_f = MODE_ROBUST
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call compute_family_scaling(&
            n_genes = n_genes,&
            n_families = n_families,&
            distances = distances,&
            gene_to_fam = gene_to_fam,&
            dscale = dscale,&
            loess_x = loess_x,&
            loess_y = loess_y,&
            indices_used = indices_used,&
            span = span,&
            degree = degree,&
            mode = mode_mode_f,&
            n_iters = n_iters,&
            low_sd_cutoff = low_sd_cutoff,&
            excluded_low_sd = excluded_low_sd,&
            ierr = ierr&
        )
    end subroutine compute_family_scaling_c

    !> summary: C-wrapper for [[tox_get_outliers(module):compute_family_scaling_expert(subroutine)]]
    !| Uses LOESS on the median/stddev of intra-family distances for scaling, regardless of orthologs.
    subroutine compute_family_scaling_expert_c(&
            n_genes,&
            n_families,&
            distances,&
            gene_to_fam,&
            dscale,&
            loess_x,&
            loess_y,&
            indices_used,&
            tmp_perm,&
            tmp_stack_left,&
            tmp_stack_right,&
            tmp_int_workspace,&
            int_workspace_size,&
            tmp_real_workspace,&
            real_workspace_size,&
            tmp_diagl,&
            tmp_weights,&
            tmp_eval_points,&
            tmp_robust_weights,&
            tmp_combined_weights,&
            tmp_residuals,&
            tmp_permutation_indices,&
            tmp_fitted_values,&
            span,&
            degree,&
            mode,&
            n_iters,&
            low_sd_cutoff,&
            excluded_low_sd,&
            tmp_means_aux,&
            ierr&
        ) bind(C, name="compute_family_scaling_expert_c")
        use tox_get_outliers, only: compute_family_scaling_expert
        use tox_loess_impl, only: MODE_PLAIN, MODE_ROBUST

        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes
        integer(c_int), intent(in), target :: n_families
            !! Total number of gene families
        integer(c_int), intent(in), target :: int_workspace_size
            !! Length of integer workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `int_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_families  |
            !! | save_factorization    | .false.     |
        integer(c_int), intent(in), target :: real_workspace_size
            !! Length of real workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `real_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_families  |
            !! | save_factorization    | .false.     |
        real(c_double), dimension(n_genes), intent(in), target :: distances
            !! Array of Euclidean distances for each gene
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Mapping of each gene to its family (1-based)
        real(c_double), dimension(n_families), intent(out), target :: dscale
            !! Array of scaling factors per family (output)
        real(c_double), dimension(n_families), intent(out), target :: loess_x
            !! Reference x-coordinates for LOESS smoothing
        real(c_double), dimension(n_families), intent(out), target :: loess_y
            !! Reference y-coordinates for LOESS smoothing
        integer(c_int), dimension(n_families), intent(out), target :: indices_used
            !! Indices of reference points used for smoothing
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_perm
            !! Permutation array for sorting gene distances
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_stack_left
            !! Stack array for left indices during sorting
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_stack_right
            !! Stack array for right indices during sorting
        integer(c_int), dimension(int_workspace_size), intent(out), target :: tmp_int_workspace
            !! Integer workspace array
        real(c_double), dimension(real_workspace_size), intent(out), target :: tmp_real_workspace
            !! Real workspace array
        real(c_double), dimension(n_families), intent(inout), target :: tmp_diagl
            !! Diagonal elements of the weight matrix
        real(c_double), dimension(n_families), intent(out), target :: tmp_weights
            !! Initial weights for LOESS
        real(c_double), dimension(n_families, 1), intent(inout), target :: tmp_eval_points
            !! Z matrix for LOESS fitting
        real(c_double), dimension(n_families), intent(out), target :: tmp_robust_weights
            !! Residuals for robust LOESS fitting
        real(c_double), dimension(n_families), intent(out), target :: tmp_combined_weights
            !! Working weights array
        real(c_double), dimension(n_families), intent(out), target :: tmp_residuals
            !! Residuals array
        integer(c_int), dimension(n_families), intent(out), target :: tmp_permutation_indices
            !! Permutation indices for robust LOESS fitting
        real(c_double), dimension(n_families), intent(out), target :: tmp_fitted_values
            !! Output array for LOESS predictions
        real(c_double), intent(in), target :: span
            !! Span parameter for LOESS smoothing, passed straight to
            !! [[tox_loess_impl(module):loess_fit_plain_impl(subroutine)]], so it is held to that
            !! procedure's own range rather than to the NaN tolerance the distance data carries.
            !! The default value is `0.7_real64`.
            !! The minimum valid value is `EPS_LOESS`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), intent(in), target :: degree
            !! Degree of the LOESS polynomial
            !! The default value is `2_int32`.
        character(len=1, kind=c_char), dimension(6), intent(in), target :: mode
            !! Mode for LOESS fitting
            !! The default value is `'robust'`.
            !!
            !! | Mode                 | Value                                            |
            !! |----------------------|--------------------------------------------------|
            !! | Plain LOESS fitting  | [[tox_loess_impl(module):MODE_PLAIN(variable)]]  |
            !! | Robust LOESS fitting | [[tox_loess_impl(module):MODE_ROBUST(variable)]] |
        integer(c_int), intent(in), target :: n_iters
            !! Number of iterations for robust LOESS fitting
            !! The default value is `3_int32`.
        real(c_double), intent(out), target :: low_sd_cutoff
            !! cutoff used to filter families with low std
        integer(c_int), dimension(n_families), intent(out), target :: excluded_low_sd
            !! Mask to save those families that have low sd
        real(c_double), dimension(n_families), intent(out), target :: tmp_means_aux
            !! Work array for saving raw means
        integer(c_int), intent(out), target :: ierr
            !! Error code
        integer(int32) :: mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(int_workspace_size)
        M_CHECK_NON_NULL(real_workspace_size)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_NON_NULL(n_iters)
        M_CHECK_NON_NULL(low_sd_cutoff)
        M_CHECK_ARRAY_NON_NULL(distances, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(dscale, n_families)
        M_CHECK_ARRAY_NON_NULL(loess_x, n_families)
        M_CHECK_ARRAY_NON_NULL(loess_y, n_families)
        M_CHECK_ARRAY_NON_NULL(indices_used, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_perm, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_stack_left, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_stack_right, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_int_workspace, int_workspace_size)
        M_CHECK_ARRAY_NON_NULL(tmp_real_workspace, real_workspace_size)
        M_CHECK_ARRAY_NON_NULL(tmp_diagl, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_eval_points, n_families * 1)
        M_CHECK_ARRAY_NON_NULL(tmp_robust_weights, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_combined_weights, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_residuals, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_indices, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_fitted_values, n_families)
        M_CHECK_ARRAY_NON_NULL(mode, 6)
        M_CHECK_ARRAY_NON_NULL(excluded_low_sd, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_means_aux, n_families)

        block
            character(len=:), allocatable :: mode_f
            call c_char_1d_as_string(mode, mode_f, ierr)
            if (is_err(ierr)) return

            select case (mode_f)
                case ("plain")
                    mode_mode_f = MODE_PLAIN
                case ("robust")
                    mode_mode_f = MODE_ROBUST
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call compute_family_scaling_expert(&
            n_genes = n_genes,&
            n_families = n_families,&
            distances = distances,&
            gene_to_fam = gene_to_fam,&
            dscale = dscale,&
            loess_x = loess_x,&
            loess_y = loess_y,&
            indices_used = indices_used,&
            tmp_perm = tmp_perm,&
            tmp_stack_left = tmp_stack_left,&
            tmp_stack_right = tmp_stack_right,&
            tmp_int_workspace = tmp_int_workspace,&
            int_workspace_size = int_workspace_size,&
            tmp_real_workspace = tmp_real_workspace,&
            real_workspace_size = real_workspace_size,&
            tmp_diagl = tmp_diagl,&
            tmp_weights = tmp_weights,&
            tmp_eval_points = tmp_eval_points,&
            tmp_robust_weights = tmp_robust_weights,&
            tmp_combined_weights = tmp_combined_weights,&
            tmp_residuals = tmp_residuals,&
            tmp_permutation_indices = tmp_permutation_indices,&
            tmp_fitted_values = tmp_fitted_values,&
            span = span,&
            degree = degree,&
            mode = mode_mode_f,&
            n_iters = n_iters,&
            low_sd_cutoff = low_sd_cutoff,&
            excluded_low_sd = excluded_low_sd,&
            tmp_means_aux = tmp_means_aux,&
            ierr = ierr&
        )
    end subroutine compute_family_scaling_expert_c

    !> summary: C-wrapper for [[tox_get_outliers(module):compute_rdi(subroutine)]]
    !| RDI = Euclidean distance / family scaling factor
    subroutine compute_rdi_c(&
            n_genes,&
            distances,&
            gene_to_fam,&
            dscale,&
            n_dscale_elements,&
            rdi,&
            sorted_rdi,&
            perm,&
            ierr&
        ) bind(C, name="compute_rdi_c")
        use tox_get_outliers, only: compute_rdi

        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes
        integer(c_int), intent(in), target :: n_dscale_elements
            !! number of elements in `dscale`
        real(c_double), dimension(n_genes), intent(in), target :: distances
            !! Array of Euclidean distances for each gene to its centroid
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Gene-to-family mapping (1-based indexing)
        real(c_double), dimension(n_dscale_elements), intent(in), target :: dscale
            !! Array of scaling factors for each family
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_genes), intent(out), target :: rdi
            !! Output array of RDI values for each gene
        real(c_double), dimension(n_genes), intent(out), target :: sorted_rdi
            !! Work array for sorting (dimension n_genes)
        integer(c_int), dimension(n_genes), intent(out), target :: perm
            !! Permutation array for sorting (dimension n_genes, should be pre-initialized with 1:n_genes)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_dscale_elements)
        M_CHECK_ARRAY_NON_NULL(distances, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(dscale, n_dscale_elements)
        M_CHECK_ARRAY_NON_NULL(rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(sorted_rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(perm, n_genes)

        call compute_rdi(&
            n_genes = n_genes,&
            distances = distances,&
            gene_to_fam = gene_to_fam,&
            dscale = dscale,&
            rdi = rdi,&
            sorted_rdi = sorted_rdi,&
            perm = perm,&
            ierr = ierr&
        )
    end subroutine compute_rdi_c

    !> summary: C-wrapper for [[tox_get_outliers(module):compute_rdi_expert(subroutine)]]
    !| RDI = Euclidean distance / family scaling factor
    subroutine compute_rdi_expert_c(&
            n_genes,&
            distances,&
            gene_to_fam,&
            dscale,&
            n_dscale_elements,&
            rdi,&
            sorted_rdi,&
            perm,&
            tmp_stack_left,&
            tmp_stack_right,&
            ierr&
        ) bind(C, name="compute_rdi_expert_c")
        use tox_get_outliers, only: compute_rdi_expert

        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes
        integer(c_int), intent(in), target :: n_dscale_elements
            !! number of elements in `dscale`
        real(c_double), dimension(n_genes), intent(in), target :: distances
            !! Array of Euclidean distances for each gene to its centroid
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Gene-to-family mapping (1-based indexing)
        real(c_double), dimension(n_dscale_elements), intent(in), target :: dscale
            !! Array of scaling factors for each family
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_genes), intent(out), target :: rdi
            !! Output array of RDI values for each gene
        real(c_double), dimension(n_genes), intent(out), target :: sorted_rdi
            !! Work array for sorting (dimension n_genes)
        integer(c_int), dimension(n_genes), intent(out), target :: perm
            !! Permutation array for sorting (dimension n_genes, should be pre-initialized with 1:n_genes)
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_stack_left
            !! Stack array for sorting (dimension n_genes)
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_stack_right
            !! Stack array for sorting (dimension n_genes)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_dscale_elements)
        M_CHECK_ARRAY_NON_NULL(distances, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(dscale, n_dscale_elements)
        M_CHECK_ARRAY_NON_NULL(rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(sorted_rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(perm, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_stack_left, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_stack_right, n_genes)

        call compute_rdi_expert(&
            n_genes = n_genes,&
            distances = distances,&
            gene_to_fam = gene_to_fam,&
            dscale = dscale,&
            rdi = rdi,&
            sorted_rdi = sorted_rdi,&
            perm = perm,&
            tmp_stack_left = tmp_stack_left,&
            tmp_stack_right = tmp_stack_right,&
            ierr = ierr&
        )
    end subroutine compute_rdi_expert_c

    !> summary: C-wrapper for [[tox_get_outliers(module):identify_outliers(subroutine)]]
    !| Expects sorted_rdi to be filtered (no negative values) and perm should be sorted in ascending order before calling.
    !| If sorted_rdi contains negatives or perm is not sorted, tmp_results may be invalid.
    subroutine identify_outliers_c(&
            n_genes,&
            rdi,&
            sorted_rdi,&
            perm,&
            is_outlier,&
            threshold,&
            quantile,&
            percentile,&
            ierr&
        ) bind(C, name="identify_outliers_c")
        use tox_get_outliers, only: identify_outliers

        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes
        real(c_double), dimension(n_genes), intent(in), target :: rdi
            !! Array of RDI values for each gene
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_genes), intent(in), target :: sorted_rdi
            !! Sorted RDI array (must be filtered to remove negatives and sorted in ascending order before calling)
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(c_int), dimension(n_genes), intent(in), target :: perm
            !! Permutation array with sorted indices
        logical(c_bool), dimension(n_genes), intent(out), target :: is_outlier
            !! Output boolean array indicating outliers
        real(c_double), intent(out), target :: threshold
            !! Output threshold value used for detection
        real(c_double), dimension(n_genes), intent(out), target :: quantile
            !! Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
            !! observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            !! Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            !! upper-tail quantile is used.
        real(c_double), intent(in), target :: percentile
            !! Percentile threshold as a fraction in [0,1] (top 5% for the default).
            !! The default value is `0.95_real64`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(threshold)
        M_CHECK_NON_NULL(percentile)
        M_CHECK_ARRAY_NON_NULL(rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(sorted_rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(perm, n_genes)
        M_CHECK_ARRAY_NON_NULL(is_outlier, n_genes)
        M_CHECK_ARRAY_NON_NULL(quantile, n_genes)

        call identify_outliers(&
            n_genes = n_genes,&
            rdi = rdi,&
            sorted_rdi = sorted_rdi,&
            perm = perm,&
            is_outlier = is_outlier,&
            threshold = threshold,&
            quantile = quantile,&
            percentile = percentile,&
            ierr = ierr&
        )
    end subroutine identify_outliers_c

    !> summary: C-wrapper for [[tox_get_outliers(module):detect_outliers(subroutine)]]
    !| Orchestrates the full pipeline: per-family scaling via
    !| [[tox_get_outliers_impl(module):compute_family_scaling_impl(subroutine)]], the RDI per gene via
    !| [[tox_get_outliers_impl(module):compute_rdi_impl(subroutine)]], then flags outliers via
    !| [[tox_get_outliers_impl(module):identify_outliers_impl(subroutine)]].
    subroutine detect_outliers_c(&
            n_genes,&
            n_families,&
            distances,&
            gene_to_fam,&
            is_outlier,&
            loess_x,&
            loess_y,&
            loess_n,&
            quantile,&
            ierr,&
            percentile&
        ) bind(C, name="detect_outliers_c")
        use tox_get_outliers, only: detect_outliers

        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes
        integer(c_int), intent(in), target :: n_families
            !! Total number of gene families
        real(c_double), dimension(n_genes), intent(in), target :: distances
            !! Array of Euclidean distances for each gene to its centroid
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Gene-to-family mapping (1-based indexing)
        logical(c_bool), dimension(n_genes), intent(out), target :: is_outlier
            !! Output boolean array indicating outliers
        real(c_double), dimension(n_families), intent(out), target :: loess_x
            !! Reference x-coordinates.
        real(c_double), dimension(n_families), intent(out), target :: loess_y
            !! Reference y-coordinates (length n_total).
        integer(c_int), dimension(n_families), intent(out), target :: loess_n
            !! Indices of reference points used for smoothing.
        real(c_double), dimension(n_genes), intent(out), target :: quantile
            !! Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
            !! observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            !! Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            !! upper-tail quantile is used.
        integer(c_int), intent(out), target :: ierr
            !! Error code
        real(c_double), intent(in), target :: percentile
            !! Percentile threshold as a fraction in [0,1] for outlier detection.
            !! The default value is `0.95_real64`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(percentile)
        M_CHECK_ARRAY_NON_NULL(distances, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(is_outlier, n_genes)
        M_CHECK_ARRAY_NON_NULL(loess_x, n_families)
        M_CHECK_ARRAY_NON_NULL(loess_y, n_families)
        M_CHECK_ARRAY_NON_NULL(loess_n, n_families)
        M_CHECK_ARRAY_NON_NULL(quantile, n_genes)

        call detect_outliers(&
            n_genes = n_genes,&
            n_families = n_families,&
            distances = distances,&
            gene_to_fam = gene_to_fam,&
            is_outlier = is_outlier,&
            loess_x = loess_x,&
            loess_y = loess_y,&
            loess_n = loess_n,&
            quantile = quantile,&
            ierr = ierr,&
            percentile = percentile&
        )
    end subroutine detect_outliers_c

    !> summary: C-wrapper for [[tox_get_outliers(module):detect_outliers_expert(subroutine)]]
    !| Orchestrates the full pipeline: per-family scaling via
    !| [[tox_get_outliers_impl(module):compute_family_scaling_impl(subroutine)]], the RDI per gene via
    !| [[tox_get_outliers_impl(module):compute_rdi_impl(subroutine)]], then flags outliers via
    !| [[tox_get_outliers_impl(module):identify_outliers_impl(subroutine)]].
    subroutine detect_outliers_expert_c(&
            n_genes,&
            n_families,&
            distances,&
            gene_to_fam,&
            tmp_perm,&
            tmp_stack_left,&
            tmp_stack_right,&
            tmp_int_workspace,&
            int_workspace_size,&
            tmp_real_workspace,&
            real_workspace_size,&
            tmp_diagl,&
            tmp_weights,&
            tmp_eval_points,&
            tmp_robust_weights,&
            tmp_combined_weights,&
            tmp_residuals,&
            tmp_permutation_indices,&
            tmp_fitted_values,&
            tmp_means_aux,&
            tmp_dscale,&
            tmp_excluded_low_sd,&
            tmp_low_sd_cutoff,&
            tmp_rdi,&
            tmp_sorted_rdi,&
            tmp_threshold,&
            is_outlier,&
            loess_x,&
            loess_y,&
            loess_n,&
            quantile,&
            ierr,&
            percentile&
        ) bind(C, name="detect_outliers_expert_c")
        use tox_get_outliers, only: detect_outliers_expert

        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes
        integer(c_int), intent(in), target :: n_families
            !! Total number of gene families
        integer(c_int), intent(in), target :: int_workspace_size
            !! Length of integer workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `int_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_families  |
            !! | save_factorization    | .false.     |
        integer(c_int), intent(in), target :: real_workspace_size
            !! Length of real workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `real_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_families  |
            !! | save_factorization    | .false.     |
        real(c_double), dimension(n_genes), intent(in), target :: distances
            !! Array of Euclidean distances for each gene to its centroid
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Gene-to-family mapping (1-based indexing)
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_perm
            !! Permutation array for sorting
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_stack_left
            !! Stack array for left indices during sorting
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_stack_right
            !! Stack array for right indices during sorting
        integer(c_int), dimension(int_workspace_size), intent(out), target :: tmp_int_workspace
            !! Integer workspace array
        real(c_double), dimension(real_workspace_size), intent(out), target :: tmp_real_workspace
            !! Real workspace array
        real(c_double), dimension(n_families), intent(out), target :: tmp_diagl
            !! Diagonal elements of the weight matrix
        real(c_double), dimension(n_families), intent(out), target :: tmp_weights
            !! Initial weights for LOESS
        real(c_double), dimension(n_families, 1), intent(out), target :: tmp_eval_points
            !! Z matrix for LOESS fitting
        real(c_double), dimension(n_families), intent(out), target :: tmp_robust_weights
            !! Residuals for robust LOESS fitting
        real(c_double), dimension(n_families), intent(out), target :: tmp_combined_weights
            !! Working weights array
        real(c_double), dimension(n_families), intent(out), target :: tmp_residuals
            !! Residuals array
        integer(c_int), dimension(n_families), intent(out), target :: tmp_permutation_indices
            !! Permutation indices for robust LOESS fitting
        real(c_double), dimension(n_families), intent(out), target :: tmp_fitted_values
            !! Output array for LOESS predictions
        real(c_double), dimension(n_families), intent(out), target :: tmp_means_aux
            !! Work array for saving raw means
        real(c_double), dimension(n_families), intent(out), target :: tmp_dscale
            !! Per-family scaling factors (intermediate, consumed by the RDI step)
        integer(c_int), dimension(n_families), intent(out), target :: tmp_excluded_low_sd
            !! Low-sd family mask (intermediate, discarded)
        real(c_double), intent(out), target :: tmp_low_sd_cutoff
            !! Low-sd cutoff (intermediate, discarded)
        real(c_double), dimension(n_genes), intent(out), target :: tmp_rdi
            !! RDI per gene (intermediate, consumed by the outlier step)
        real(c_double), dimension(n_genes), intent(out), target :: tmp_sorted_rdi
            !! Sorted RDI (intermediate, consumed by the outlier step)
        real(c_double), intent(out), target :: tmp_threshold
            !! Detection threshold (intermediate, discarded)
        logical(c_bool), dimension(n_genes), intent(out), target :: is_outlier
            !! Output boolean array indicating outliers
        real(c_double), dimension(n_families), intent(out), target :: loess_x
            !! Reference x-coordinates.
        real(c_double), dimension(n_families), intent(out), target :: loess_y
            !! Reference y-coordinates (length n_total).
        integer(c_int), dimension(n_families), intent(out), target :: loess_n
            !! Indices of reference points used for smoothing.
        real(c_double), dimension(n_genes), intent(out), target :: quantile
            !! Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
            !! observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            !! Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            !! upper-tail quantile is used.
        integer(c_int), intent(out), target :: ierr
            !! Error code
        real(c_double), intent(in), target :: percentile
            !! Percentile threshold as a fraction in [0,1] for outlier detection.
            !! The default value is `0.95_real64`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(int_workspace_size)
        M_CHECK_NON_NULL(real_workspace_size)
        M_CHECK_NON_NULL(tmp_low_sd_cutoff)
        M_CHECK_NON_NULL(tmp_threshold)
        M_CHECK_NON_NULL(percentile)
        M_CHECK_ARRAY_NON_NULL(distances, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_perm, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_stack_left, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_stack_right, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_int_workspace, int_workspace_size)
        M_CHECK_ARRAY_NON_NULL(tmp_real_workspace, real_workspace_size)
        M_CHECK_ARRAY_NON_NULL(tmp_diagl, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_eval_points, n_families * 1)
        M_CHECK_ARRAY_NON_NULL(tmp_robust_weights, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_combined_weights, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_residuals, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_indices, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_fitted_values, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_means_aux, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_dscale, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_excluded_low_sd, n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_sorted_rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(is_outlier, n_genes)
        M_CHECK_ARRAY_NON_NULL(loess_x, n_families)
        M_CHECK_ARRAY_NON_NULL(loess_y, n_families)
        M_CHECK_ARRAY_NON_NULL(loess_n, n_families)
        M_CHECK_ARRAY_NON_NULL(quantile, n_genes)

        call detect_outliers_expert(&
            n_genes = n_genes,&
            n_families = n_families,&
            distances = distances,&
            gene_to_fam = gene_to_fam,&
            tmp_perm = tmp_perm,&
            tmp_stack_left = tmp_stack_left,&
            tmp_stack_right = tmp_stack_right,&
            tmp_int_workspace = tmp_int_workspace,&
            int_workspace_size = int_workspace_size,&
            tmp_real_workspace = tmp_real_workspace,&
            real_workspace_size = real_workspace_size,&
            tmp_diagl = tmp_diagl,&
            tmp_weights = tmp_weights,&
            tmp_eval_points = tmp_eval_points,&
            tmp_robust_weights = tmp_robust_weights,&
            tmp_combined_weights = tmp_combined_weights,&
            tmp_residuals = tmp_residuals,&
            tmp_permutation_indices = tmp_permutation_indices,&
            tmp_fitted_values = tmp_fitted_values,&
            tmp_means_aux = tmp_means_aux,&
            tmp_dscale = tmp_dscale,&
            tmp_excluded_low_sd = tmp_excluded_low_sd,&
            tmp_low_sd_cutoff = tmp_low_sd_cutoff,&
            tmp_rdi = tmp_rdi,&
            tmp_sorted_rdi = tmp_sorted_rdi,&
            tmp_threshold = tmp_threshold,&
            is_outlier = is_outlier,&
            loess_x = loess_x,&
            loess_y = loess_y,&
            loess_n = loess_n,&
            quantile = quantile,&
            ierr = ierr,&
            percentile = percentile&
        )
    end subroutine detect_outliers_expert_c

end module tox_get_outliers_c
#endif
