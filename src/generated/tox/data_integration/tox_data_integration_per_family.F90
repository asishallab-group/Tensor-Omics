#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation per family
!|
!| The JSD value for a sub-neighborhood -- typically the genes of one family -- obtained by
!| driving the same pipeline over a masked set of neighbors. Answers whether two studies are
!| compatible *for this family*, which the global figure can hide either way.
!|
!| Generated from [[tox_data_integration_per_family_impl(module)]]; do not edit -- regenerate instead.
module tox_data_integration_per_family
    use tox_data_integration_per_family_impl, only: fjct_compute_contribution_scores_impl, fjct_compute_jsd_impl, fjct_compute_masked_jsd_impl
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size, validate_in_range_int
    use tox_errors, only: validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: fjct_compute_jsd
    public :: fjct_compute_jsd_expert
    public :: fjct_compute_masked_jsd
    public :: fjct_compute_masked_jsd_expert
    public :: fjct_compute_contribution_scores

contains

    !> summary: Validates its inputs, prepares what [[tox_data_integration_per_family_impl(module):fjct_compute_jsd_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_per_family(module):fjct_compute_jsd_expert]] to prepare it yourself.
    !| Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts
    !| the residual samples to the genes belonging to the family `family_idx`.
    pure subroutine fjct_compute_jsd(&
            family_idx,&
            gene_to_family_S1,&
            gene_to_family_S2,&
            n_genes_S1,&
            n_genes_S2,&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            neighborhood_genes_S1,&
            neighborhood_genes_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            n_bins,&
            shared_residual_range,&
            js_divergences,&
            included_n_reps_S1,&
            included_n_reps_S2,&
            total_included_n_reps,&
            global_js_divergence,&
            weights,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes_S1
            !! Number of genes in study 1
        integer(int32), intent(in) :: n_genes_S2
            !! Number of genes in study 2
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: family_idx
            !! Index of the family that should be analyzed
            !! The minimum valid value is `1_int32`.
        integer(int32), dimension(n_genes_S1), intent(in) :: gene_to_family_S1
            !! Mapping for study 1: Each index (gene) holds the index of its family
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_genes_S2), intent(in) :: gene_to_family_S2
            !! Mapping for study 2: Each index (gene) holds the index of its family
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S1
            !! Indices of the selected neighborhood genes of study 1
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes_S1`.
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S2
            !! Indices of the selected neighborhood genes of study 2
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes_S2`.
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
        integer(int32), intent(out) :: total_included_n_reps
            !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical(c_bool), dimension(:, :), allocatable :: tmp_neighbor_mask_S1
        logical(c_bool), dimension(:, :), allocatable :: tmp_neighbor_mask_S2
        real(real64), dimension(:, :), allocatable :: tmp_pmf_S1
        real(real64), dimension(:, :), allocatable :: tmp_pmf_S2
        integer(int32), dimension(:, :), allocatable :: tmp_counts

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(family_idx, ierr, arg_pos=1_int32, min=1_int32)
        call validate_dimension_size(n_genes_S1, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_genes_S2, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_reps_S1, ierr, arg_pos=10_int32)
        call validate_dimension_size(n_reps_S2, ierr, arg_pos=11_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=12_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=13_int32)
        call validate_in_range_real(shared_residual_range, ierr, arg_pos=15_int32, min=0.0_real64)
        call validate_all_in_range_int(gene_to_family_S1, n_genes_S1, ierr, arg_pos=2_int32, min=0_int32)
        call validate_all_in_range_int(gene_to_family_S2, n_genes_S2, ierr, arg_pos=3_int32, min=0_int32)
        call validate_all_in_range_real(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points, ierr, arg_pos=6_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points, ierr, arg_pos=7_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(neighborhood_genes_S1, n_neighbors * n_points, ierr, arg_pos=8_int32, min=1_int32, max=n_genes_S1)
        call validate_all_in_range_int(neighborhood_genes_S2, n_neighbors * n_points, ierr, arg_pos=9_int32, min=1_int32, max=n_genes_S2)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_neighbor_mask_S1(n_neighbors, n_points))
        M_ALLOCATE(tmp_neighbor_mask_S2(n_neighbors, n_points))
        M_ALLOCATE(tmp_pmf_S1(n_points, n_bins))
        M_ALLOCATE(tmp_pmf_S2(n_points, n_bins))
        M_ALLOCATE(tmp_counts(n_points, n_bins))

        call fjct_compute_jsd_impl(&
            family_idx = family_idx,&
            gene_to_family_S1 = gene_to_family_S1,&
            gene_to_family_S2 = gene_to_family_S2,&
            n_genes_S1 = n_genes_S1,&
            n_genes_S2 = n_genes_S2,&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            neighborhood_genes_S1 = neighborhood_genes_S1,&
            neighborhood_genes_S2 = neighborhood_genes_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            js_divergences = js_divergences,&
            included_n_reps_S1 = included_n_reps_S1,&
            included_n_reps_S2 = included_n_reps_S2,&
            total_included_n_reps = total_included_n_reps,&
            global_js_divergence = global_js_divergence,&
            weights = weights,&
            tmp_neighbor_mask_S1 = tmp_neighbor_mask_S1,&
            tmp_neighbor_mask_S2 = tmp_neighbor_mask_S2,&
            tmp_pmf_S1 = tmp_pmf_S1,&
            tmp_pmf_S2 = tmp_pmf_S2,&
            tmp_counts = tmp_counts&
        )
    end subroutine fjct_compute_jsd

    !> summary: Validates its inputs, then calls [[tox_data_integration_per_family_impl(module):fjct_compute_jsd_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_per_family(module):fjct_compute_jsd]] does both.
    !| Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts
    !| the residual samples to the genes belonging to the family `family_idx`.
    pure subroutine fjct_compute_jsd_expert(&
            family_idx,&
            gene_to_family_S1,&
            gene_to_family_S2,&
            n_genes_S1,&
            n_genes_S2,&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            neighborhood_genes_S1,&
            neighborhood_genes_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            n_bins,&
            shared_residual_range,&
            js_divergences,&
            included_n_reps_S1,&
            included_n_reps_S2,&
            total_included_n_reps,&
            global_js_divergence,&
            weights,&
            tmp_neighbor_mask_S1,&
            tmp_neighbor_mask_S2,&
            tmp_pmf_S1,&
            tmp_pmf_S2,&
            tmp_counts,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes_S1
            !! Number of genes in study 1
        integer(int32), intent(in) :: n_genes_S2
            !! Number of genes in study 2
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies
        integer(int32), intent(in) :: family_idx
            !! Index of the family that should be analyzed
            !! The minimum valid value is `1_int32`.
        integer(int32), dimension(n_genes_S1), intent(in) :: gene_to_family_S1
            !! Mapping for study 1: Each index (gene) holds the index of its family
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_genes_S2), intent(in) :: gene_to_family_S2
            !! Mapping for study 2: Each index (gene) holds the index of its family
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S1
            !! Indices of the selected neighborhood genes of study 1
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes_S1`.
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S2
            !! Indices of the selected neighborhood genes of study 2
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes_S2`.
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
        integer(int32), intent(out) :: total_included_n_reps
            !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        logical(c_bool), dimension(n_neighbors, n_points), intent(out) :: tmp_neighbor_mask_S1
            !! Work array for the mask selecting study 1's neighbors that belong to `family_idx`
        logical(c_bool), dimension(n_neighbors, n_points), intent(out) :: tmp_neighbor_mask_S2
            !! Work array for the mask selecting study 2's neighbors that belong to `family_idx`
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S1
            !! Work array for study 1's normalized histogram counts
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S2
            !! Work array for study 2's normalized histogram counts
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Work array for the histogram counts
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(family_idx, ierr, arg_pos=1_int32, min=1_int32)
        call validate_dimension_size(n_genes_S1, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_genes_S2, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_reps_S1, ierr, arg_pos=10_int32)
        call validate_dimension_size(n_reps_S2, ierr, arg_pos=11_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=12_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=13_int32)
        call validate_dimension_size(n_bins, ierr, arg_pos=14_int32)
        call validate_in_range_real(shared_residual_range, ierr, arg_pos=15_int32, min=0.0_real64)
        call validate_all_in_range_int(gene_to_family_S1, n_genes_S1, ierr, arg_pos=2_int32, min=0_int32)
        call validate_all_in_range_int(gene_to_family_S2, n_genes_S2, ierr, arg_pos=3_int32, min=0_int32)
        call validate_all_in_range_real(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points, ierr, arg_pos=6_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points, ierr, arg_pos=7_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(neighborhood_genes_S1, n_neighbors * n_points, ierr, arg_pos=8_int32, min=1_int32, max=n_genes_S1)
        call validate_all_in_range_int(neighborhood_genes_S2, n_neighbors * n_points, ierr, arg_pos=9_int32, min=1_int32, max=n_genes_S2)
        if (is_err(ierr)) return
#endif

        call fjct_compute_jsd_impl(&
            family_idx = family_idx,&
            gene_to_family_S1 = gene_to_family_S1,&
            gene_to_family_S2 = gene_to_family_S2,&
            n_genes_S1 = n_genes_S1,&
            n_genes_S2 = n_genes_S2,&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            neighborhood_genes_S1 = neighborhood_genes_S1,&
            neighborhood_genes_S2 = neighborhood_genes_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            js_divergences = js_divergences,&
            included_n_reps_S1 = included_n_reps_S1,&
            included_n_reps_S2 = included_n_reps_S2,&
            total_included_n_reps = total_included_n_reps,&
            global_js_divergence = global_js_divergence,&
            weights = weights,&
            tmp_neighbor_mask_S1 = tmp_neighbor_mask_S1,&
            tmp_neighbor_mask_S2 = tmp_neighbor_mask_S2,&
            tmp_pmf_S1 = tmp_pmf_S1,&
            tmp_pmf_S2 = tmp_pmf_S2,&
            tmp_counts = tmp_counts&
        )
    end subroutine fjct_compute_jsd_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_per_family_impl(module):fjct_compute_masked_jsd_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_per_family(module):fjct_compute_masked_jsd_expert]] to prepare it yourself.
    !| Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts the
    !| residual samples to the neighbors selected by `neighbor_mask_S1`/`neighbor_mask_S2`. Typically
    !| those are all neighbors belonging to one gene family, which is what `fjct_compute_jsd` builds
    !| the masks for from a family index.
    pure subroutine fjct_compute_masked_jsd(&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            neighbor_mask_S1,&
            neighbor_mask_S2,&
            n_bins,&
            shared_residual_range,&
            js_divergences,&
            included_n_reps_S1,&
            included_n_reps_S2,&
            total_included_n_reps,&
            global_js_divergence,&
            weights,&
            pmf_S1,&
            pmf_S2,&
            ierr&
        )
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        logical(c_bool), dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S1
            !! Mask selecting the neighbors of study 1 to include
        logical(c_bool), dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S2
            !! Mask selecting the neighbors of study 2 to include
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
        integer(int32), intent(out) :: total_included_n_reps
            !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S1
            !! Normalized histogram counts for study 1
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S2
            !! Normalized histogram counts for study 2
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:, :), allocatable :: tmp_counts

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_reps_S1, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_reps_S2, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=6_int32)
        call validate_dimension_size(n_bins, ierr, arg_pos=9_int32)
        call validate_in_range_real(shared_residual_range, ierr, arg_pos=10_int32, min=0.0_real64)
        call validate_all_in_range_real(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points, ierr, arg_pos=2_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_counts(n_points, n_bins))

        call fjct_compute_masked_jsd_impl(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            neighbor_mask_S1 = neighbor_mask_S1,&
            neighbor_mask_S2 = neighbor_mask_S2,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            js_divergences = js_divergences,&
            included_n_reps_S1 = included_n_reps_S1,&
            included_n_reps_S2 = included_n_reps_S2,&
            total_included_n_reps = total_included_n_reps,&
            global_js_divergence = global_js_divergence,&
            weights = weights,&
            pmf_S1 = pmf_S1,&
            pmf_S2 = pmf_S2,&
            tmp_counts = tmp_counts&
        )
    end subroutine fjct_compute_masked_jsd

    !> summary: Validates its inputs, then calls [[tox_data_integration_per_family_impl(module):fjct_compute_masked_jsd_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_per_family(module):fjct_compute_masked_jsd]] does both.
    !| Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts the
    !| residual samples to the neighbors selected by `neighbor_mask_S1`/`neighbor_mask_S2`. Typically
    !| those are all neighbors belonging to one gene family, which is what `fjct_compute_jsd` builds
    !| the masks for from a family index.
    pure subroutine fjct_compute_masked_jsd_expert(&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            neighbor_mask_S1,&
            neighbor_mask_S2,&
            n_bins,&
            shared_residual_range,&
            js_divergences,&
            included_n_reps_S1,&
            included_n_reps_S2,&
            total_included_n_reps,&
            global_js_divergence,&
            weights,&
            pmf_S1,&
            pmf_S2,&
            tmp_counts,&
            ierr&
        )
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        logical(c_bool), dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S1
            !! Mask selecting the neighbors of study 1 to include
        logical(c_bool), dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S2
            !! Mask selecting the neighbors of study 2 to include
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
        integer(int32), intent(out) :: total_included_n_reps
            !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S1
            !! Normalized histogram counts for study 1
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S2
            !! Normalized histogram counts for study 2
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Work array for the histogram counts
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_reps_S1, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_reps_S2, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=6_int32)
        call validate_dimension_size(n_bins, ierr, arg_pos=9_int32)
        call validate_in_range_real(shared_residual_range, ierr, arg_pos=10_int32, min=0.0_real64)
        call validate_all_in_range_real(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points, ierr, arg_pos=2_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        call fjct_compute_masked_jsd_impl(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            neighbor_mask_S1 = neighbor_mask_S1,&
            neighbor_mask_S2 = neighbor_mask_S2,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            js_divergences = js_divergences,&
            included_n_reps_S1 = included_n_reps_S1,&
            included_n_reps_S2 = included_n_reps_S2,&
            total_included_n_reps = total_included_n_reps,&
            global_js_divergence = global_js_divergence,&
            weights = weights,&
            pmf_S1 = pmf_S1,&
            pmf_S2 = pmf_S2,&
            tmp_counts = tmp_counts&
        )
    end subroutine fjct_compute_masked_jsd_expert

    !> summary: Validates its inputs, then calls [[tox_data_integration_per_family_impl(module):fjct_compute_contribution_scores_impl]].
    !| Combines
    !|
    !| 1. how divergent the family is between the studies (`global_js_divergences`), and
    !| 2. how much residual support the family has overall (`total_included_n_reps_per_f`),
    !|
    !| using the outputs of `fjct_compute_jsd`, collected for the analyzed sub-neighborhoods.
    pure subroutine fjct_compute_contribution_scores(&
            global_js_divergences,&
            total_included_n_reps_per_f,&
            k_families,&
            support_weights,&
            contribution_scores,&
            ierr&
        )
        integer(int32), intent(in) :: k_families
            !! Number of sub-neighborhoods analyzed
        real(real64), dimension(k_families), intent(in) :: global_js_divergences
            !! Per-sub-neighborhood weighted global JSD
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(k_families), intent(in) :: total_included_n_reps_per_f
            !! Per-sub-neighborhood `total_included_n_reps`
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(k_families), intent(out) :: support_weights
            !! Per-sub-neighborhood calculated support weight (ratio between its `total_included_n_reps` and `sum(total_included_n_reps_per_f)`, zero if there were no replicates included at all)
        real(real64), dimension(k_families), intent(out) :: contribution_scores
            !! Per-sub-neighborhood calculated contribution ( \( support\_weights_i * global\_js\_divergences_i \) )
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(k_families, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(global_js_divergences, k_families, ierr, arg_pos=1_int32, min=0.0_real64)
        call validate_all_in_range_int(total_included_n_reps_per_f, k_families, ierr, arg_pos=2_int32, min=0_int32)
        if (is_err(ierr)) return
#endif

        call fjct_compute_contribution_scores_impl(&
            global_js_divergences = global_js_divergences,&
            total_included_n_reps_per_f = total_included_n_reps_per_f,&
            k_families = k_families,&
            support_weights = support_weights,&
            contribution_scores = contribution_scores&
        )
    end subroutine fjct_compute_contribution_scores

end module tox_data_integration_per_family
