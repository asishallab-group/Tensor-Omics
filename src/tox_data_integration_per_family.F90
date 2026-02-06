#include "macros.h"

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation per family
!|
!| This module implements subroutines to obtain the JSD value from neighborhood residuals obtained from [[tox_data_integration_preprocessing(submodule)]]
!| for specific sub-neighborhoods by using the pipeline from [[tox_data_integration_jsd(submodule)]].
submodule (tox_data_integration) tox_data_integration_per_family
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, validate_dimension_size, validate_in_range_int, validate_all_in_range_int, validate_all_in_range_real, is_err, ERR_ALLOC_FAIL, set_err, validate_in_range_real
    use f42_utils, only: is_close
    implicit none

contains

    !> Computes the family-level compatibility score `global_js_divergence` between two studies for a single gene family (`family_idx`), by reusing the same conditioning-on-mean-expression pipeline as the global gJCT, but restricting residual samples to genes belonging to the specified family
    pure module subroutine fjct_compute_jsd_alloc(family_idx, gene_to_family_S1, gene_to_family_S2, n_genes_S1, n_genes_S2, neighborhood_residuals_S1, neighborhood_residuals_S2, &
            neighborhood_genes_S1, neighborhood_genes_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, js_divergences, &
            included_n_reps_S1, included_n_reps_S2, total_included_n_reps, global_js_divergence, weights, ierr &
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
        integer(int32), dimension(n_genes_S1), intent(in) :: gene_to_family_S1
            !! Mapping for study 1: Each index (gene) holds the index of its family
        integer(int32), dimension(n_genes_S2), intent(in) :: gene_to_family_S2
            !! Mapping for study 2: Each index (gene) holds the index of its family
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S1
            !! Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_data_integration(module):construct_neighborhoods(interface)]]
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S2
            !! Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_data_integration(module):construct_neighborhoods(interface)]]
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(int32), intent(out) :: total_included_n_reps
            !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_point, i_neighbor
        logical, dimension(:, :), allocatable :: neighbor_mask_S1, neighbor_mask_S2
        real(real64), dimension(:, :), allocatable :: pmf_S1, pmf_S2
        integer(int32), dimension(:, :), allocatable :: tmp_counts

        call set_ok(ierr)

        call validate_dimension_size(n_neighbors, ierr)
        call validate_in_range_int(family_idx, ierr, min=1_int32)
        call validate_all_in_range_int(gene_to_family_S1, n_genes_S1, ierr, min=0_int32)
        call validate_all_in_range_int(gene_to_family_S2, n_genes_S2, ierr, min=0_int32)
        call validate_all_in_range_int(neighborhood_genes_S1, size(neighborhood_genes_S1, kind=int32), ierr, min=1_int32, max=n_genes_S1)
        call validate_all_in_range_int(neighborhood_genes_S2, size(neighborhood_genes_S2, kind=int32), ierr, min=1_int32, max=n_genes_S2)

        if (is_err(ierr)) return

        M_ALLOCATE(neighbor_mask_S1(n_neighbors, n_points))
        M_ALLOCATE(neighbor_mask_S2(n_neighbors, n_points))
        M_ALLOCATE(pmf_S1(n_points, n_bins))
        M_ALLOCATE(pmf_S2(n_points, n_bins))
        M_ALLOCATE(tmp_counts(n_points, n_bins))

        ! Set up mask for filtered analysis -> only include neighbors being part of the family
        do concurrent (i_point = 1:n_points)
            do concurrent (i_neighbor = 1:n_neighbors) shared(family_idx, neighbor_mask_S1, neighbor_mask_S2, gene_to_family_S1, neighborhood_genes_S1, gene_to_family_S2, neighborhood_genes_S2)
                neighbor_mask_S1(i_neighbor, i_point) = gene_to_family_S1(neighborhood_genes_S1(i_neighbor, i_point)) == family_idx
                neighbor_mask_S2(i_neighbor, i_point) = gene_to_family_S2(neighborhood_genes_S2(i_neighbor, i_point)) == family_idx
            end do
        end do

        call fjct_compute_jsd(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, neighbor_mask_S1, neighbor_mask_S2, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, total_included_n_reps, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, ierr)
    end subroutine fjct_compute_jsd_alloc

    !> Computes the compatibility score `global_js_divergence` between two studies per sub-neighborhood/family for a single gene family (`family_idx`), by reusing the same conditioning-on-mean-expression pipeline as the global gJCT, but restricting residual samples to genes belonging to the specified family
    pure module subroutine fjct_compute_jsd(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, neighbor_mask_S1, neighbor_mask_S2, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, total_included_n_reps, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, ierr)
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        logical, dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical, dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(int32), intent(out) :: total_included_n_reps
            !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S1
            !! Absolute counts of a residual per bin (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S2
            !! Absolute counts of a residual per bin (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_reps_S1, ierr)
        call validate_dimension_size(n_reps_S2, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_in_range_int(n_bins, ierr, min=1_int32)
        call validate_in_range_real(shared_residual_range, ierr, min=0.0_real64)

        if (is_err(ierr)) return

        call jct_compute_jsd_pipeline_helper(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, neighbor_mask_S1, neighbor_mask_S2)
        total_included_n_reps = sum(included_n_reps_S1) + sum(included_n_reps_S2)
    end subroutine fjct_compute_jsd

    !> Computes the per-family/per-sub-neighborhood contribution score that combines
    !|
    !| 1. how divergent the family is between the studies (``), and
    !| 2. how much residual support the family has overall (),
    !|
    !| using the outputs from [[tox_data_integration_per_family(module):fjct_compute_jsd(subroutine)]], collected for the analyzed sub-neighborhoods.
    pure module subroutine fjct_compute_contribution_scores(global_js_divergences, total_included_n_reps_per_f, k_families, support_weights, contribution_scores, ierr)
        integer(int32), intent(in) :: k_families
            !! Number of sub-neighborhoods analyzed
        integer(int32), dimension(k_families), intent(in) :: total_included_n_reps_per_f
            !! Per-sub-neighborhood `total_included_n_reps`
        real(real64), dimension(k_families), intent(in) :: global_js_divergences
            !! Per-sub-neighborhood weighted global JSD
        real(real64), dimension(k_families), intent(out) :: support_weights
            !! Per-sub-neighborhood calculated support weight (ratio between its `total_included_n_reps` and `sum(total_included_n_reps_per_f)`, zero if there were no replicates included at all)
        real(real64), dimension(k_families), intent(out) :: contribution_scores
            !! Per-sub-neighborhood calculated contribution ( \( support\_weights_i * global\_js\_divergences_i \) )
        integer(int32), intent(out), target :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(k_families, ierr)
        call validate_all_in_range_int(total_included_n_reps_per_f, k_families, ierr, min=0_int32)
        call validate_all_in_range_real(global_js_divergences, k_families, ierr, min=0.0_real64)
        call validate_all_in_range_real(support_weights, k_families, ierr, min=0.0_real64)
        call validate_all_in_range_real(contribution_scores, k_families, ierr, min=0.0_real64)

        if (is_err(ierr)) return

        call fjct_compute_contribution_scores_helper(global_js_divergences, total_included_n_reps_per_f, k_families, support_weights, contribution_scores)
    end subroutine fjct_compute_contribution_scores

    !> (no input validation) Computes the per-family/per-sub-neighborhood contribution score that combines
    !|
    !| 1. how divergent the family is between the studies (``), and
    !| 2. how much residual support the family has overall (),
    !|
    !| using the outputs from [[tox_data_integration_per_family(module):fjct_compute_jsd(subroutine)]], collected for the analyzed sub-neighborhoods.
    pure module subroutine fjct_compute_contribution_scores_helper(global_js_divergences, total_included_n_reps_per_f, k_families, support_weights, contribution_scores)
        integer(int32), intent(in) :: k_families
            !! Number of sub-neighborhoods analyzed
        integer(int32), dimension(k_families), intent(in) :: total_included_n_reps_per_f
            !! Per-sub-neighborhood `total_included_n_reps`
        real(real64), dimension(k_families), intent(in) :: global_js_divergences
            !! Per-sub-neighborhood weighted global JSD
        real(real64), dimension(k_families), intent(out) :: support_weights
            !! Per-sub-neighborhood calculated support weight (ratio between its `total_included_n_reps` and `sum(total_included_n_reps_per_f)`, zero if there were no replicates included at all)
        real(real64), dimension(k_families), intent(out) :: contribution_scores
            !! Per-sub-neighborhood calculated contribution ( \( support\_weights_i * global\_js\_divergences_i \) )

        integer(int32) :: i_family
        real(real64) :: total_included_n_reps

        total_included_n_reps = real(sum(total_included_n_reps_per_f), kind=real64)

        if (is_close(total_included_n_reps, 0.0_real64)) then
            support_weights = 0.0_real64
            contribution_scores = 0.0_real64
        end if

        do concurrent (i_family = 1:k_families) shared(support_weights, total_included_n_reps_per_f, total_included_n_reps, contribution_scores, global_js_divergences)
            support_weights(i_family) = real(total_included_n_reps_per_f(i_family), kind=real64) / total_included_n_reps

            contribution_scores(i_family) = support_weights(i_family) * global_js_divergences(i_family)
        end do
    end subroutine fjct_compute_contribution_scores_helper

end submodule tox_data_integration_per_family

!> C-compatible wrapper for [[tox_data_integration(module):fjct_compute_jsd_alloc(interface)]]
pure subroutine fjct_compute_jsd_c( &
        family_idx, gene_to_family_S1, gene_to_family_S2, n_genes_S1, n_genes_S2, &
        neighborhood_residuals_S1, neighborhood_residuals_S2, &
        neighborhood_genes_S1, neighborhood_genes_S2, &
        n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, &
        js_divergences, included_n_reps_S1, included_n_reps_S2, total_included_n_reps, &
        global_js_divergence, weights, ierr &
    ) bind(C, name="fjct_compute_jsd_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_data_integration, only: fjct_compute_jsd_alloc
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes_S1
        !! Number of genes in study 1
    integer(c_int), intent(in), target :: n_genes_S2
        !! Number of genes in study 2
    integer(c_int), intent(in), target :: n_reps_S1
        !! Number of replicates in study 1
    integer(c_int), intent(in), target :: n_reps_S2
        !! Number of replicates in study 2
    integer(c_int), intent(in), target :: n_neighbors
        !! Number of neighbors in the studies
    integer(c_int), intent(in), target :: n_points
        !! Number of reference points in the studies
    integer(c_int), intent(in), target :: family_idx
        !! Index of the family that should be analyzed

    integer(c_int), dimension(n_genes_S1), intent(in), target :: gene_to_family_S1
        !! Mapping for study 1: Each index (gene) holds the index of its family
    integer(c_int), dimension(n_genes_S2), intent(in), target :: gene_to_family_S2
        !! Mapping for study 2: Each index (gene) holds the index of its family

    real(c_double), dimension(n_reps_S1,n_neighbors,n_points), intent(in), target :: neighborhood_residuals_S1
        !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]),
        !! NaN is explicitly allowed for missing values
    real(c_double), dimension(n_reps_S2,n_neighbors,n_points), intent(in), target :: neighborhood_residuals_S2
        !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]),
        !! NaN is explicitly allowed for missing values

    integer(c_int), dimension(n_neighbors,n_points), intent(in), target :: neighborhood_genes_S1
        !! Indices of selected neighborhood genes, obtained from `neighborhood_indices`
        !! of [[tox_data_integration(module):construct_neighborhoods(interface)]]
    integer(c_int), dimension(n_neighbors,n_points), intent(in), target :: neighborhood_genes_S2
        !! Indices of selected neighborhood genes, obtained from `neighborhood_indices`
        !! of [[tox_data_integration(module):construct_neighborhoods(interface)]]

    integer(c_int), intent(in), target :: n_bins
        !! Number of equally sized histogram bins used for the studies in
        !! [[tox_data_integration(module):build_residual_histograms(interface)]]

    real(c_double), intent(in), target :: shared_residual_range
        !! Computed residual range for both studies, from
        !! [[tox_data_integration(module):determine_shared_residual_range(interface)]]

    real(c_double), dimension(n_points), intent(out), target :: js_divergences
        !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
    integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S1
        !! Count of non-NaN residuals (included ones) in study 1
    integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S2
        !! Count of non-NaN residuals (included ones) in study 2
    integer(c_int), intent(out), target :: total_included_n_reps
        !! Total number of included replicates from both studies
    real(c_double), intent(out), target :: global_js_divergence
        !! Weighted global Jensen-Shannon divergence
    real(c_double), dimension(n_points), intent(out), target :: weights
        !! Weights used for calculating the global weighted Jensen-Shannon divergence

    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes_S1)
    M_CHECK_NON_NULL(n_genes_S2)
    M_CHECK_NON_NULL(n_reps_S1)
    M_CHECK_NON_NULL(n_reps_S2)
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(family_idx)
    M_CHECK_NON_NULL(gene_to_family_S1)
    M_CHECK_NON_NULL(gene_to_family_S2)
    M_CHECK_NON_NULL(neighborhood_residuals_S1)
    M_CHECK_NON_NULL(neighborhood_residuals_S2)
    M_CHECK_NON_NULL(neighborhood_genes_S1)
    M_CHECK_NON_NULL(neighborhood_genes_S2)
    M_CHECK_NON_NULL(n_bins)
    M_CHECK_NON_NULL(shared_residual_range)
    M_CHECK_NON_NULL(js_divergences)
    M_CHECK_NON_NULL(included_n_reps_S1)
    M_CHECK_NON_NULL(included_n_reps_S2)
    M_CHECK_NON_NULL(total_included_n_reps)
    M_CHECK_NON_NULL(global_js_divergence)
    M_CHECK_NON_NULL(weights)

    call fjct_compute_jsd_alloc( &
        family_idx, gene_to_family_S1, gene_to_family_S2, n_genes_S1, n_genes_S2, &
        neighborhood_residuals_S1, neighborhood_residuals_S2, &
        neighborhood_genes_S1, neighborhood_genes_S2, &
        n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, &
        js_divergences, included_n_reps_S1, included_n_reps_S2, total_included_n_reps, &
        global_js_divergence, weights, ierr )

end subroutine fjct_compute_jsd_c

!> C-compatible wrapper for [[tox_data_integration(module):fjct_compute_jsd(interface)]]
pure subroutine fjct_compute_jsd_expert_c( &
        neighborhood_residuals_S1, neighborhood_residuals_S2, &
        n_reps_S1, n_reps_S2, n_neighbors, n_points, &
        neighbor_mask_S1, neighbor_mask_S2, &
        n_bins, shared_residual_range, &
        js_divergences, included_n_reps_S1, included_n_reps_S2, &
        total_included_n_reps, global_js_divergence, weights, &
        pmf_S1, pmf_S2, tmp_counts, ierr &
    ) bind(C, name="fjct_compute_jsd_expert_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_data_integration, only: fjct_compute_jsd
    use tox_conversions, only: c_int_as_logical
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_reps_S1
        !! Number of replicates in study 1
    integer(c_int), intent(in), target :: n_reps_S2
        !! Number of replicates in study 2
    integer(c_int), intent(in), target :: n_neighbors
        !! Number of neighbors in the studies
    integer(c_int), intent(in), target :: n_points
        !! Number of reference points in the studies

    real(c_double), dimension(n_reps_S1,n_neighbors,n_points), intent(in), target :: neighborhood_residuals_S1
        !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]),
        !! NaN is explicitly allowed for missing values
    real(c_double), dimension(n_reps_S2,n_neighbors,n_points), intent(in), target :: neighborhood_residuals_S2
        !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]),
        !! NaN is explicitly allowed for missing values

    integer(c_int), dimension(n_neighbors,n_points), intent(in), target :: neighbor_mask_S1
        !! Optional mask to exclude specific neighbors from study 1
    integer(c_int), dimension(n_neighbors,n_points), intent(in), target :: neighbor_mask_S2
        !! Optional mask to exclude specific neighbors from study 2

    integer(c_int), intent(in), target :: n_bins
        !! Number of equally sized histogram bins
    real(c_double), intent(in), target :: shared_residual_range
        !! Shared residual range

    real(c_double), dimension(n_points), intent(out), target :: js_divergences
        !! Jensen-Shannon divergence per reference point
    integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S1
        !! Included replicates S1
    integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S2
        !! Included replicates S2
    integer(c_int), intent(out), target :: total_included_n_reps
        !! Total included replicates
    real(c_double), intent(out), target :: global_js_divergence
        !! Weighted global JSD
    real(c_double), dimension(n_points), intent(out), target :: weights
        !! Weights

    real(c_double), dimension(n_points,n_bins), intent(out), target :: pmf_S1
        !! PMF S1
    real(c_double), dimension(n_points,n_bins), intent(out), target :: pmf_S2
        !! PMF S2
    integer(c_int), dimension(n_points,n_bins), intent(out), target :: tmp_counts
        !! Temporary histogram counts

    integer(c_int), intent(out), target :: ierr
        !! Error code

    logical, dimension(n_neighbors,n_points) :: mask1_f, mask2_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_reps_S1)
    M_CHECK_NON_NULL(n_reps_S2)
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(neighborhood_residuals_S1)
    M_CHECK_NON_NULL(neighborhood_residuals_S2)
    M_CHECK_NON_NULL(neighbor_mask_S1)
    M_CHECK_NON_NULL(neighbor_mask_S2)
    M_CHECK_NON_NULL(n_bins)
    M_CHECK_NON_NULL(shared_residual_range)
    M_CHECK_NON_NULL(js_divergences)
    M_CHECK_NON_NULL(included_n_reps_S1)
    M_CHECK_NON_NULL(included_n_reps_S2)
    M_CHECK_NON_NULL(total_included_n_reps)
    M_CHECK_NON_NULL(global_js_divergence)
    M_CHECK_NON_NULL(weights)
    M_CHECK_NON_NULL(pmf_S1)
    M_CHECK_NON_NULL(pmf_S2)
    M_CHECK_NON_NULL(tmp_counts)

    call c_int_as_logical(neighbor_mask_S1, mask1_f)
    call c_int_as_logical(neighbor_mask_S2, mask2_f)

    call fjct_compute_jsd( &
        neighborhood_residuals_S1, neighborhood_residuals_S2, &
        n_reps_S1, n_reps_S2, n_neighbors, n_points, &
        mask1_f, mask2_f, &
        n_bins, shared_residual_range, &
        js_divergences, included_n_reps_S1, included_n_reps_S2, &
        total_included_n_reps, global_js_divergence, weights, &
        pmf_S1, pmf_S2, tmp_counts, ierr )

end subroutine fjct_compute_jsd_expert_c

!> C-compatible wrapper for [[tox_data_integration(module):fjct_compute_contribution_scores(interface)]]
pure subroutine fjct_compute_contribution_scores_c( &
        global_js_divergences, total_included_n_reps_per_f, k_families, &
        support_weights, contribution_scores, ierr &
    ) bind(C, name="fjct_compute_contribution_scores_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_data_integration, only: fjct_compute_contribution_scores
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: k_families
        !! Number of sub-neighborhoods analyzed
    integer(c_int), dimension(k_families), intent(in), target :: total_included_n_reps_per_f
        !! Per-sub-neighborhood `total_included_n_reps`
    real(c_double), dimension(k_families), intent(in), target :: global_js_divergences
        !! Per-sub-neighborhood weighted global JSD
    real(c_double), dimension(k_families), intent(out), target :: support_weights
        !! Per-sub-neighborhood calculated support weight (ratio between its
        !! `total_included_n_reps` and `sum(total_included_n_reps_per_f)`,
        !! zero if there were no replicates included at all)
    real(c_double), dimension(k_families), intent(out), target :: contribution_scores
        !! Per-sub-neighborhood calculated contribution
        !! ( \( support\_weights_i * global\_js\_divergences_i \) )
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(k_families)
    M_CHECK_NON_NULL(total_included_n_reps_per_f)
    M_CHECK_NON_NULL(global_js_divergences)
    M_CHECK_NON_NULL(support_weights)
    M_CHECK_NON_NULL(contribution_scores)

    call fjct_compute_contribution_scores( &
        global_js_divergences, total_included_n_reps_per_f, k_families, &
        support_weights, contribution_scores, ierr )
end subroutine fjct_compute_contribution_scores_c
