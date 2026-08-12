#include <src/macros.h>

!> Gene outliers, from how far each gene sits from its family's centroid.
!|
!| The pipeline is three steps, each callable on its own. `compute_rdi` turns raw distances into
!| a relative distance index, scaled per family so families of different spread are comparable.
!| `compute_family_scaling` fits that scaling with LOESS against family size. `identify_outliers`
!| applies the threshold and reports which genes exceed it.
!|
!| `detect_outliers` runs all three in one call, and is the entry point to reach for first.
!|
!| Generated from [[tox_get_outliers_impl(module)]]; do not edit -- regenerate instead.
module tox_get_outliers
    use tox_get_outliers_impl, only: compute_family_scaling_impl, compute_rdi_impl, detect_outliers_impl, identify_outliers_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_loess_impl, only: EPS_LOESS, MODE_PLAIN, MODE_ROBUST, tox_loess_required_workspace
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    use tox_errors, only: clear_err_arg_pos, set_err, set_err_once, validate_all_in_range_real
    use tox_errors, only: validate_dimension_size, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: compute_family_scaling
    public :: compute_family_scaling_expert
    public :: compute_rdi
    public :: compute_rdi_expert
    public :: identify_outliers
    public :: detect_outliers
    public :: detect_outliers_expert

contains

    !> summary: Validates its inputs, prepares what [[tox_get_outliers_impl(module):compute_family_scaling_impl]] needs, then calls it. The entry point to reach for first; see [[tox_get_outliers(module):compute_family_scaling_expert]] to prepare it yourself.
    !| Uses LOESS on the median/stddev of intra-family distances for scaling, regardless of orthologs.
    subroutine compute_family_scaling(&
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
        )
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        real(real64), dimension(n_genes), intent(in) :: distances
            !! Array of Euclidean distances for each gene
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Mapping of each gene to its family (1-based)
        real(real64), dimension(n_families), intent(out) :: dscale
            !! Array of scaling factors per family (output)
        real(real64), dimension(n_families), intent(out) :: loess_x
            !! Reference x-coordinates for LOESS smoothing
        real(real64), dimension(n_families), intent(out) :: loess_y
            !! Reference y-coordinates for LOESS smoothing
        integer(int32), dimension(n_families), intent(out) :: indices_used
            !! Indices of reference points used for smoothing
        real(real64), intent(in), optional :: span
            !! Span parameter for LOESS smoothing, passed straight to
            !! [[tox_loess_impl(module):loess_fit_plain_impl(subroutine)]], so it is held to that
            !! procedure's own range rather than to the NaN tolerance the distance data carries.
            !! The default value is `0.7_real64`.
            !! The minimum valid value is `EPS_LOESS`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), intent(in), optional :: degree
            !! Degree of the LOESS polynomial
            !! The default value is `2_int32`.
        integer(int32), intent(in), optional :: mode
            !! Mode for LOESS fitting
            !! The default value is `1_int32`.
            !!
            !! | Mode                 | Value                                            |
            !! |----------------------|--------------------------------------------------|
            !! | Plain LOESS fitting  | [[tox_loess_impl(module):MODE_PLAIN(variable)]]  |
            !! | Robust LOESS fitting | [[tox_loess_impl(module):MODE_ROBUST(variable)]] |
        integer(int32), intent(in), optional :: n_iters
            !! Number of iterations for robust LOESS fitting
            !! The default value is `3_int32`.
        real(real64), intent(out) :: low_sd_cutoff
            !! cutoff used to filter families with low std
        integer(int32), dimension(n_families), intent(out) :: excluded_low_sd
            !! Mask to save those families that have low sd
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), dimension(:), allocatable :: tmp_perm
        integer(int32), dimension(:), allocatable :: tmp_stack_left
        integer(int32), dimension(:), allocatable :: tmp_stack_right
        integer(int32), dimension(:), allocatable :: tmp_int_workspace
        integer(int32) :: int_workspace_size
        real(real64), dimension(:), allocatable :: tmp_real_workspace
        integer(int32) :: real_workspace_size
        real(real64), dimension(:), allocatable :: tmp_diagl
        real(real64), dimension(:), allocatable :: tmp_weights
        real(real64), dimension(:, :), allocatable :: tmp_eval_points
        real(real64), dimension(:), allocatable :: tmp_robust_weights
        real(real64), dimension(:), allocatable :: tmp_combined_weights
        real(real64), dimension(:), allocatable :: tmp_residuals
        integer(int32), dimension(:), allocatable :: tmp_permutation_indices
        real(real64), dimension(:), allocatable :: tmp_fitted_values
        real(real64), dimension(:), allocatable :: tmp_means_aux

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=2_int32)
        call validate_in_range_real(span, ierr, arg_pos=9_int32, min=EPS_LOESS, max=1.0_real64)
        if (present(mode)) then; if (mode /= MODE_PLAIN .and. mode /= MODE_ROBUST) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=11_int32); end if
        if (is_err(ierr)) return
#endif

        call tox_loess_required_workspace(&
            n_dim = 1_int32,&
            max_neighborhood_size = n_families,&
            int_workspace_size = int_workspace_size,&
            real_workspace_size = real_workspace_size,&
            save_factorization = .false.&
        )
        M_ALLOCATE(tmp_perm(n_genes))
        M_ALLOCATE(tmp_stack_left(n_genes))
        M_ALLOCATE(tmp_stack_right(n_genes))
        M_ALLOCATE(tmp_int_workspace(int_workspace_size))
        M_ALLOCATE(tmp_real_workspace(real_workspace_size))
        M_ALLOCATE(tmp_diagl(n_families))
        M_ALLOCATE(tmp_weights(n_families))
        M_ALLOCATE(tmp_eval_points(n_families, 1))
        M_ALLOCATE(tmp_robust_weights(n_families))
        M_ALLOCATE(tmp_combined_weights(n_families))
        M_ALLOCATE(tmp_residuals(n_families))
        M_ALLOCATE(tmp_permutation_indices(n_families))
        M_ALLOCATE(tmp_fitted_values(n_families))
        M_ALLOCATE(tmp_means_aux(n_families))

        call compute_family_scaling_impl(&
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
            mode = mode,&
            n_iters = n_iters,&
            low_sd_cutoff = low_sd_cutoff,&
            excluded_low_sd = excluded_low_sd,&
            tmp_means_aux = tmp_means_aux,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine compute_family_scaling

    !> summary: Validates its inputs, then calls [[tox_get_outliers_impl(module):compute_family_scaling_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_get_outliers(module):compute_family_scaling]] does both.
    !| Uses LOESS on the median/stddev of intra-family distances for scaling, regardless of orthologs.
    subroutine compute_family_scaling_expert(&
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
        )
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        integer(int32), intent(in) :: int_workspace_size
            !! Length of integer workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `int_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_families  |
            !! | save_factorization    | .false.     |
        integer(int32), intent(in) :: real_workspace_size
            !! Length of real workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `real_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_families  |
            !! | save_factorization    | .false.     |
        real(real64), dimension(n_genes), intent(in) :: distances
            !! Array of Euclidean distances for each gene
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Mapping of each gene to its family (1-based)
        real(real64), dimension(n_families), intent(out) :: dscale
            !! Array of scaling factors per family (output)
        real(real64), dimension(n_families), intent(out) :: loess_x
            !! Reference x-coordinates for LOESS smoothing
        real(real64), dimension(n_families), intent(out) :: loess_y
            !! Reference y-coordinates for LOESS smoothing
        integer(int32), dimension(n_families), intent(out) :: indices_used
            !! Indices of reference points used for smoothing
        integer(int32), dimension(n_genes), intent(out) :: tmp_perm
            !! Permutation array for sorting gene distances
        integer(int32), dimension(n_genes), intent(out) :: tmp_stack_left
            !! Stack array for left indices during sorting
        integer(int32), dimension(n_genes), intent(out) :: tmp_stack_right
            !! Stack array for right indices during sorting
        integer(int32), dimension(int_workspace_size), intent(out) :: tmp_int_workspace
            !! Integer workspace array
        real(real64), dimension(real_workspace_size), intent(out) :: tmp_real_workspace
            !! Real workspace array
        real(real64), dimension(n_families), intent(inout) :: tmp_diagl
            !! Diagonal elements of the weight matrix
        real(real64), dimension(n_families), intent(out) :: tmp_weights
            !! Initial weights for LOESS
        real(real64), dimension(n_families, 1), intent(inout) :: tmp_eval_points
            !! Z matrix for LOESS fitting
        real(real64), dimension(n_families), intent(out) :: tmp_robust_weights
            !! Residuals for robust LOESS fitting
        real(real64), dimension(n_families), intent(out) :: tmp_combined_weights
            !! Working weights array
        real(real64), dimension(n_families), intent(out) :: tmp_residuals
            !! Residuals array
        integer(int32), dimension(n_families), intent(out) :: tmp_permutation_indices
            !! Permutation indices for robust LOESS fitting
        real(real64), dimension(n_families), intent(out) :: tmp_fitted_values
            !! Output array for LOESS predictions
        real(real64), intent(in), optional :: span
            !! Span parameter for LOESS smoothing, passed straight to
            !! [[tox_loess_impl(module):loess_fit_plain_impl(subroutine)]], so it is held to that
            !! procedure's own range rather than to the NaN tolerance the distance data carries.
            !! The default value is `0.7_real64`.
            !! The minimum valid value is `EPS_LOESS`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), intent(in), optional :: degree
            !! Degree of the LOESS polynomial
            !! The default value is `2_int32`.
        integer(int32), intent(in), optional :: mode
            !! Mode for LOESS fitting
            !! The default value is `1_int32`.
            !!
            !! | Mode                 | Value                                            |
            !! |----------------------|--------------------------------------------------|
            !! | Plain LOESS fitting  | [[tox_loess_impl(module):MODE_PLAIN(variable)]]  |
            !! | Robust LOESS fitting | [[tox_loess_impl(module):MODE_ROBUST(variable)]] |
        integer(int32), intent(in), optional :: n_iters
            !! Number of iterations for robust LOESS fitting
            !! The default value is `3_int32`.
        real(real64), intent(out) :: low_sd_cutoff
            !! cutoff used to filter families with low std
        integer(int32), dimension(n_families), intent(out) :: excluded_low_sd
            !! Mask to save those families that have low sd
        real(real64), dimension(n_families), intent(out) :: tmp_means_aux
            !! Work array for saving raw means
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=2_int32)
        call validate_dimension_size(int_workspace_size, ierr, arg_pos=13_int32)
        call validate_dimension_size(real_workspace_size, ierr, arg_pos=15_int32)
        call validate_in_range_real(span, ierr, arg_pos=24_int32, min=EPS_LOESS, max=1.0_real64)
        call validate_all_in_range_real(tmp_diagl, n_families, ierr, arg_pos=16_int32)
        call validate_all_in_range_real(tmp_eval_points, n_families * 1, ierr, arg_pos=18_int32)
        if (present(mode)) then; if (mode /= MODE_PLAIN .and. mode /= MODE_ROBUST) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=26_int32); end if
        if (is_err(ierr)) return
#endif

        call compute_family_scaling_impl(&
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
            mode = mode,&
            n_iters = n_iters,&
            low_sd_cutoff = low_sd_cutoff,&
            excluded_low_sd = excluded_low_sd,&
            tmp_means_aux = tmp_means_aux,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine compute_family_scaling_expert

    !> summary: Validates its inputs, prepares what [[tox_get_outliers_impl(module):compute_rdi_impl]] needs, then calls it. The entry point to reach for first; see [[tox_get_outliers(module):compute_rdi_expert]] to prepare it yourself.
    !| RDI = Euclidean distance / family scaling factor
    pure subroutine compute_rdi(&
            n_genes,&
            distances,&
            gene_to_fam,&
            dscale,&
            rdi,&
            sorted_rdi,&
            perm,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        real(real64), dimension(n_genes), intent(in) :: distances
            !! Array of Euclidean distances for each gene to its centroid
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene-to-family mapping (1-based indexing)
        real(real64), dimension(:), intent(in) :: dscale
            !! Array of scaling factors for each family
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_genes), intent(out) :: rdi
            !! Output array of RDI values for each gene
        real(real64), dimension(n_genes), intent(out) :: sorted_rdi
            !! Work array for sorting (dimension n_genes)
        integer(int32), dimension(n_genes), intent(out) :: perm
            !! Permutation array for sorting (dimension n_genes, should be pre-initialized with 1:n_genes)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: tmp_stack_left
        integer(int32), dimension(:), allocatable :: tmp_stack_right

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_stack_left(n_genes))
        M_ALLOCATE(tmp_stack_right(n_genes))

        call compute_rdi_impl(&
            n_genes = n_genes,&
            distances = distances,&
            gene_to_fam = gene_to_fam,&
            dscale = dscale,&
            rdi = rdi,&
            sorted_rdi = sorted_rdi,&
            perm = perm,&
            tmp_stack_left = tmp_stack_left,&
            tmp_stack_right = tmp_stack_right&
        )
    end subroutine compute_rdi

    !> summary: Validates its inputs, then calls [[tox_get_outliers_impl(module):compute_rdi_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_get_outliers(module):compute_rdi]] does both.
    !| RDI = Euclidean distance / family scaling factor
    pure subroutine compute_rdi_expert(&
            n_genes,&
            distances,&
            gene_to_fam,&
            dscale,&
            rdi,&
            sorted_rdi,&
            perm,&
            tmp_stack_left,&
            tmp_stack_right,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        real(real64), dimension(n_genes), intent(in) :: distances
            !! Array of Euclidean distances for each gene to its centroid
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene-to-family mapping (1-based indexing)
        real(real64), dimension(:), intent(in) :: dscale
            !! Array of scaling factors for each family
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_genes), intent(out) :: rdi
            !! Output array of RDI values for each gene
        real(real64), dimension(n_genes), intent(out) :: sorted_rdi
            !! Work array for sorting (dimension n_genes)
        integer(int32), dimension(n_genes), intent(out) :: perm
            !! Permutation array for sorting (dimension n_genes, should be pre-initialized with 1:n_genes)
        integer(int32), dimension(n_genes), intent(out) :: tmp_stack_left
            !! Stack array for sorting (dimension n_genes)
        integer(int32), dimension(n_genes), intent(out) :: tmp_stack_right
            !! Stack array for sorting (dimension n_genes)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        call compute_rdi_impl(&
            n_genes = n_genes,&
            distances = distances,&
            gene_to_fam = gene_to_fam,&
            dscale = dscale,&
            rdi = rdi,&
            sorted_rdi = sorted_rdi,&
            perm = perm,&
            tmp_stack_left = tmp_stack_left,&
            tmp_stack_right = tmp_stack_right&
        )
    end subroutine compute_rdi_expert

    !> summary: Validates its inputs, then calls [[tox_get_outliers_impl(module):identify_outliers_impl]].
    !| Expects sorted_rdi to be filtered (no negative values) and perm should be sorted in ascending order before calling.
    !| If sorted_rdi contains negatives or perm is not sorted, tmp_results may be invalid.
    pure subroutine identify_outliers(&
            n_genes,&
            rdi,&
            sorted_rdi,&
            perm,&
            is_outlier,&
            threshold,&
            quantile,&
            percentile,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        real(real64), dimension(n_genes), intent(in) :: rdi
            !! Array of RDI values for each gene
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_genes), intent(in) :: sorted_rdi
            !! Sorted RDI array (must be filtered to remove negatives and sorted in ascending order before calling)
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(int32), dimension(n_genes), intent(in) :: perm
            !! Permutation array with sorted indices
        logical, dimension(n_genes), intent(out) :: is_outlier
            !! Output boolean array indicating outliers
        real(real64), intent(out) :: threshold
            !! Output threshold value used for detection
        real(real64), dimension(n_genes), intent(out) :: quantile
            !! Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
            !! observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            !! Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            !! upper-tail quantile is used.
        real(real64), intent(in), optional :: percentile
            !! Percentile threshold (top 5% for the default).
            !! The default value is `95.0_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_in_range_real(percentile, ierr, arg_pos=8_int32)
        if (is_err(ierr)) return
#endif

        call identify_outliers_impl(&
            n_genes = n_genes,&
            rdi = rdi,&
            sorted_rdi = sorted_rdi,&
            perm = perm,&
            is_outlier = is_outlier,&
            threshold = threshold,&
            quantile = quantile,&
            percentile = percentile&
        )
    end subroutine identify_outliers

    !> summary: Validates its inputs, prepares what [[tox_get_outliers_impl(module):detect_outliers_impl]] needs, then calls it. The entry point to reach for first; see [[tox_get_outliers(module):detect_outliers_expert]] to prepare it yourself.
    !| Orchestrates the full pipeline: per-family scaling via
    !| [[tox_get_outliers_impl(module):compute_family_scaling_impl(subroutine)]], the RDI per gene via
    !| [[tox_get_outliers_impl(module):compute_rdi_impl(subroutine)]], then flags outliers via
    !| [[tox_get_outliers_impl(module):identify_outliers_impl(subroutine)]].
    subroutine detect_outliers(&
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
        )
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        real(real64), dimension(n_genes), intent(in) :: distances
            !! Array of Euclidean distances for each gene to its centroid
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene-to-family mapping (1-based indexing)
        logical, dimension(n_genes), intent(out) :: is_outlier
            !! Output boolean array indicating outliers
        real(real64), dimension(n_families), intent(out) :: loess_x
            !! Reference x-coordinates.
        real(real64), dimension(n_families), intent(out) :: loess_y
            !! Reference y-coordinates (length n_total).
        integer(int32), dimension(n_families), intent(out) :: loess_n
            !! Indices of reference points used for smoothing.
        real(real64), dimension(n_genes), intent(out) :: quantile
            !! Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
            !! observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            !! Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            !! upper-tail quantile is used.
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), intent(in), optional :: percentile
            !! Percentile threshold for outlier detection.
            !! The default value is `95.0_real64`.
        integer(int32), dimension(:), allocatable :: tmp_perm
        integer(int32), dimension(:), allocatable :: tmp_stack_left
        integer(int32), dimension(:), allocatable :: tmp_stack_right
        integer(int32), dimension(:), allocatable :: tmp_int_workspace
        integer(int32) :: int_workspace_size
        real(real64), dimension(:), allocatable :: tmp_real_workspace
        integer(int32) :: real_workspace_size
        real(real64), dimension(:), allocatable :: tmp_diagl
        real(real64), dimension(:), allocatable :: tmp_weights
        real(real64), dimension(:, :), allocatable :: tmp_eval_points
        real(real64), dimension(:), allocatable :: tmp_robust_weights
        real(real64), dimension(:), allocatable :: tmp_combined_weights
        real(real64), dimension(:), allocatable :: tmp_residuals
        integer(int32), dimension(:), allocatable :: tmp_permutation_indices
        real(real64), dimension(:), allocatable :: tmp_fitted_values
        real(real64), dimension(:), allocatable :: tmp_means_aux
        real(real64), dimension(:), allocatable :: tmp_dscale
        integer(int32), dimension(:), allocatable :: tmp_excluded_low_sd
        real(real64) :: tmp_low_sd_cutoff
        real(real64), dimension(:), allocatable :: tmp_rdi
        real(real64), dimension(:), allocatable :: tmp_sorted_rdi
        real(real64) :: tmp_threshold

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=2_int32)
        call validate_in_range_real(percentile, ierr, arg_pos=11_int32)
        if (is_err(ierr)) return
#endif

        call tox_loess_required_workspace(&
            n_dim = 1_int32,&
            max_neighborhood_size = n_families,&
            int_workspace_size = int_workspace_size,&
            real_workspace_size = real_workspace_size,&
            save_factorization = .false.&
        )
        M_ALLOCATE(tmp_perm(n_genes))
        M_ALLOCATE(tmp_stack_left(n_genes))
        M_ALLOCATE(tmp_stack_right(n_genes))
        M_ALLOCATE(tmp_int_workspace(int_workspace_size))
        M_ALLOCATE(tmp_real_workspace(real_workspace_size))
        M_ALLOCATE(tmp_diagl(n_families))
        M_ALLOCATE(tmp_weights(n_families))
        M_ALLOCATE(tmp_eval_points(n_families, 1))
        M_ALLOCATE(tmp_robust_weights(n_families))
        M_ALLOCATE(tmp_combined_weights(n_families))
        M_ALLOCATE(tmp_residuals(n_families))
        M_ALLOCATE(tmp_permutation_indices(n_families))
        M_ALLOCATE(tmp_fitted_values(n_families))
        M_ALLOCATE(tmp_means_aux(n_families))
        M_ALLOCATE(tmp_dscale(n_families))
        M_ALLOCATE(tmp_excluded_low_sd(n_families))
        M_ALLOCATE(tmp_rdi(n_genes))
        M_ALLOCATE(tmp_sorted_rdi(n_genes))

        call detect_outliers_impl(&
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
        call clear_err_arg_pos(ierr)
    end subroutine detect_outliers

    !> summary: Validates its inputs, then calls [[tox_get_outliers_impl(module):detect_outliers_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_get_outliers(module):detect_outliers]] does both.
    !| Orchestrates the full pipeline: per-family scaling via
    !| [[tox_get_outliers_impl(module):compute_family_scaling_impl(subroutine)]], the RDI per gene via
    !| [[tox_get_outliers_impl(module):compute_rdi_impl(subroutine)]], then flags outliers via
    !| [[tox_get_outliers_impl(module):identify_outliers_impl(subroutine)]].
    subroutine detect_outliers_expert(&
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
        )
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        integer(int32), intent(in) :: int_workspace_size
            !! Length of integer workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `int_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_families  |
            !! | save_factorization    | .false.     |
        integer(int32), intent(in) :: real_workspace_size
            !! Length of real workspace.
            !! It is *VERY IMPORTANT* to compute this argument from the `real_workspace_size` output produced by [[tox_loess_impl(module):tox_loess_required_workspace]].
            !!
            !! | Producer input        | Supplied by |
            !! |-----------------------|-------------|
            !! | n_dim                 | 1_int32     |
            !! | max_neighborhood_size | n_families  |
            !! | save_factorization    | .false.     |
        real(real64), dimension(n_genes), intent(in) :: distances
            !! Array of Euclidean distances for each gene to its centroid
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene-to-family mapping (1-based indexing)
        integer(int32), dimension(n_genes), intent(out) :: tmp_perm
            !! Permutation array for sorting
        integer(int32), dimension(n_genes), intent(out) :: tmp_stack_left
            !! Stack array for left indices during sorting
        integer(int32), dimension(n_genes), intent(out) :: tmp_stack_right
            !! Stack array for right indices during sorting
        integer(int32), dimension(int_workspace_size), intent(out) :: tmp_int_workspace
            !! Integer workspace array
        real(real64), dimension(real_workspace_size), intent(out) :: tmp_real_workspace
            !! Real workspace array
        real(real64), dimension(n_families), intent(out) :: tmp_diagl
            !! Diagonal elements of the weight matrix
        real(real64), dimension(n_families), intent(out) :: tmp_weights
            !! Initial weights for LOESS
        real(real64), dimension(n_families, 1), intent(out) :: tmp_eval_points
            !! Z matrix for LOESS fitting
        real(real64), dimension(n_families), intent(out) :: tmp_robust_weights
            !! Residuals for robust LOESS fitting
        real(real64), dimension(n_families), intent(out) :: tmp_combined_weights
            !! Working weights array
        real(real64), dimension(n_families), intent(out) :: tmp_residuals
            !! Residuals array
        integer(int32), dimension(n_families), intent(out) :: tmp_permutation_indices
            !! Permutation indices for robust LOESS fitting
        real(real64), dimension(n_families), intent(out) :: tmp_fitted_values
            !! Output array for LOESS predictions
        real(real64), dimension(n_families), intent(out) :: tmp_means_aux
            !! Work array for saving raw means
        real(real64), dimension(n_families), intent(out) :: tmp_dscale
            !! Per-family scaling factors (intermediate, consumed by the RDI step)
        integer(int32), dimension(n_families), intent(out) :: tmp_excluded_low_sd
            !! Low-sd family mask (intermediate, discarded)
        real(real64), intent(out) :: tmp_low_sd_cutoff
            !! Low-sd cutoff (intermediate, discarded)
        real(real64), dimension(n_genes), intent(out) :: tmp_rdi
            !! RDI per gene (intermediate, consumed by the outlier step)
        real(real64), dimension(n_genes), intent(out) :: tmp_sorted_rdi
            !! Sorted RDI (intermediate, consumed by the outlier step)
        real(real64), intent(out) :: tmp_threshold
            !! Detection threshold (intermediate, discarded)
        logical, dimension(n_genes), intent(out) :: is_outlier
            !! Output boolean array indicating outliers
        real(real64), dimension(n_families), intent(out) :: loess_x
            !! Reference x-coordinates.
        real(real64), dimension(n_families), intent(out) :: loess_y
            !! Reference y-coordinates (length n_total).
        integer(int32), dimension(n_families), intent(out) :: loess_n
            !! Indices of reference points used for smoothing.
        real(real64), dimension(n_genes), intent(out) :: quantile
            !! Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
            !! observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            !! Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            !! upper-tail quantile is used.
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), intent(in), optional :: percentile
            !! Percentile threshold for outlier detection.
            !! The default value is `95.0_real64`.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=2_int32)
        call validate_dimension_size(int_workspace_size, ierr, arg_pos=9_int32)
        call validate_dimension_size(real_workspace_size, ierr, arg_pos=11_int32)
        call validate_in_range_real(percentile, ierr, arg_pos=33_int32)
        if (is_err(ierr)) return
#endif

        call detect_outliers_impl(&
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
        call clear_err_arg_pos(ierr)
    end subroutine detect_outliers_expert

end module tox_get_outliers
