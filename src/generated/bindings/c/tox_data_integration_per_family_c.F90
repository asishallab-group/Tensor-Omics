#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_per_family(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation per family
!|
!| The JSD value for a sub-neighborhood -- typically the genes of one family -- obtained by
!| driving the same pipeline over a masked set of neighbors. Answers whether two studies are
!| compatible *for this family*, which the global figure can hide either way.
module tox_data_integration_per_family_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_ALLOC_FAIL
    M_IMPLICIT_NONE
    private

    public :: fjct_compute_jsd_c
    public :: fjct_compute_jsd_expert_c
    public :: fjct_compute_masked_jsd_c
    public :: fjct_compute_masked_jsd_expert_c
    public :: fjct_compute_contribution_scores_c

contains

    !> summary: C-wrapper for [[tox_data_integration_per_family(module):fjct_compute_jsd(subroutine)]]
    !| Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts
    !| the residual samples to the genes belonging to the family `family_idx`.
    subroutine fjct_compute_jsd_c(&
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
        ) bind(C, name="fjct_compute_jsd_c")
        use tox_data_integration_per_family, only: fjct_compute_jsd

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
            !! The minimum valid value is `1_int32`.
        integer(c_int), dimension(n_genes_S1), intent(in), target :: gene_to_family_S1
            !! Mapping for study 1: Each index (gene) holds the index of its family
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_genes_S2), intent(in), target :: gene_to_family_S2
            !! Mapping for study 2: Each index (gene) holds the index of its family
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_neighbors, n_points), intent(in), target :: neighborhood_genes_S1
            !! Indices of the selected neighborhood genes of study 1
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes_S1`.
        integer(c_int), dimension(n_neighbors, n_points), intent(in), target :: neighborhood_genes_S2
            !! Indices of the selected neighborhood genes of study 2
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes_S2`.
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_points), intent(out), target :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
        integer(c_int), intent(out), target :: total_included_n_reps
            !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        real(c_double), intent(out), target :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(c_double), dimension(n_points), intent(out), target :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(family_idx)
        M_CHECK_NON_NULL(n_genes_S1)
        M_CHECK_NON_NULL(n_genes_S2)
        M_CHECK_NON_NULL(n_reps_S1)
        M_CHECK_NON_NULL(n_reps_S2)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(total_included_n_reps)
        M_CHECK_NON_NULL(global_js_divergence)
        M_CHECK_ARRAY_NON_NULL(gene_to_family_S1, n_genes_S1)
        M_CHECK_ARRAY_NON_NULL(gene_to_family_S2, n_genes_S2)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_genes_S1, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_genes_S2, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S1, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S2, n_points)
        M_CHECK_ARRAY_NON_NULL(weights, n_points)

        call fjct_compute_jsd(&
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
            ierr = ierr&
        )
    end subroutine fjct_compute_jsd_c

    !> summary: C-wrapper for [[tox_data_integration_per_family(module):fjct_compute_jsd_expert(subroutine)]]
    !| Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts
    !| the residual samples to the genes belonging to the family `family_idx`.
    subroutine fjct_compute_jsd_expert_c(&
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
        ) bind(C, name="fjct_compute_jsd_expert_c")
        use tox_data_integration_per_family, only: fjct_compute_jsd_expert

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
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins used for the studies
        integer(c_int), intent(in), target :: family_idx
            !! Index of the family that should be analyzed
            !! The minimum valid value is `1_int32`.
        integer(c_int), dimension(n_genes_S1), intent(in), target :: gene_to_family_S1
            !! Mapping for study 1: Each index (gene) holds the index of its family
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_genes_S2), intent(in), target :: gene_to_family_S2
            !! Mapping for study 2: Each index (gene) holds the index of its family
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_neighbors, n_points), intent(in), target :: neighborhood_genes_S1
            !! Indices of the selected neighborhood genes of study 1
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes_S1`.
        integer(c_int), dimension(n_neighbors, n_points), intent(in), target :: neighborhood_genes_S2
            !! Indices of the selected neighborhood genes of study 2
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes_S2`.
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_points), intent(out), target :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
        integer(c_int), intent(out), target :: total_included_n_reps
            !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        real(c_double), intent(out), target :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(c_double), dimension(n_points), intent(out), target :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        logical(c_bool), dimension(n_neighbors, n_points), intent(out), target :: tmp_neighbor_mask_S1
            !! Work array for the mask selecting study 1's neighbors that belong to `family_idx`
        logical(c_bool), dimension(n_neighbors, n_points), intent(out), target :: tmp_neighbor_mask_S2
            !! Work array for the mask selecting study 2's neighbors that belong to `family_idx`
        real(c_double), dimension(n_points, n_bins), intent(out), target :: tmp_pmf_S1
            !! Work array for study 1's normalized histogram counts
        real(c_double), dimension(n_points, n_bins), intent(out), target :: tmp_pmf_S2
            !! Work array for study 2's normalized histogram counts
        integer(c_int), dimension(n_points, n_bins), intent(out), target :: tmp_counts
            !! Work array for the histogram counts
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(:, :), allocatable :: tmp_neighbor_mask_S1_f
        logical, dimension(:, :), allocatable :: tmp_neighbor_mask_S2_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(family_idx)
        M_CHECK_NON_NULL(n_genes_S1)
        M_CHECK_NON_NULL(n_genes_S2)
        M_CHECK_NON_NULL(n_reps_S1)
        M_CHECK_NON_NULL(n_reps_S2)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(total_included_n_reps)
        M_CHECK_NON_NULL(global_js_divergence)
        M_CHECK_ARRAY_NON_NULL(gene_to_family_S1, n_genes_S1)
        M_CHECK_ARRAY_NON_NULL(gene_to_family_S2, n_genes_S2)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_genes_S1, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_genes_S2, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S1, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S2, n_points)
        M_CHECK_ARRAY_NON_NULL(weights, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_neighbor_mask_S1, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_neighbor_mask_S2, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_pmf_S1, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_pmf_S2, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, n_points * n_bins)

        M_ALLOCATE(tmp_neighbor_mask_S1_f(n_neighbors, n_points))
        M_ALLOCATE(tmp_neighbor_mask_S2_f(n_neighbors, n_points))

        call fjct_compute_jsd_expert(&
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
            tmp_neighbor_mask_S1 = tmp_neighbor_mask_S1_f,&
            tmp_neighbor_mask_S2 = tmp_neighbor_mask_S2_f,&
            tmp_pmf_S1 = tmp_pmf_S1,&
            tmp_pmf_S2 = tmp_pmf_S2,&
            tmp_counts = tmp_counts,&
            ierr = ierr&
        )

        tmp_neighbor_mask_S1 = tmp_neighbor_mask_S1_f
        tmp_neighbor_mask_S2 = tmp_neighbor_mask_S2_f
    end subroutine fjct_compute_jsd_expert_c

    !> summary: C-wrapper for [[tox_data_integration_per_family(module):fjct_compute_masked_jsd(subroutine)]]
    !| Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts the
    !| residual samples to the neighbors selected by `neighbor_mask_S1`/`neighbor_mask_S2`. Typically
    !| those are all neighbors belonging to one gene family, which is what `fjct_compute_jsd` builds
    !| the masks for from a family index.
    subroutine fjct_compute_masked_jsd_c(&
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
        ) bind(C, name="fjct_compute_masked_jsd_c")
        use tox_data_integration_per_family, only: fjct_compute_masked_jsd

        integer(c_int), intent(in), target :: n_reps_S1
            !! Number of replicates in study 1
        integer(c_int), intent(in), target :: n_reps_S2
            !! Number of replicates in study 2
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors in the studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), target :: neighbor_mask_S1
            !! Mask selecting the neighbors of study 1 to include
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), target :: neighbor_mask_S2
            !! Mask selecting the neighbors of study 2 to include
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_points), intent(out), target :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
        integer(c_int), intent(out), target :: total_included_n_reps
            !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        real(c_double), intent(out), target :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(c_double), dimension(n_points), intent(out), target :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        real(c_double), dimension(n_points, n_bins), intent(out), target :: pmf_S1
            !! Normalized histogram counts for study 1
        real(c_double), dimension(n_points, n_bins), intent(out), target :: pmf_S2
            !! Normalized histogram counts for study 2
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(:, :), allocatable :: neighbor_mask_S1_f
        logical, dimension(:, :), allocatable :: neighbor_mask_S2_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_reps_S1)
        M_CHECK_NON_NULL(n_reps_S2)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(total_included_n_reps)
        M_CHECK_NON_NULL(global_js_divergence)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighbor_mask_S1, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighbor_mask_S2, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S1, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S2, n_points)
        M_CHECK_ARRAY_NON_NULL(weights, n_points)
        M_CHECK_ARRAY_NON_NULL(pmf_S1, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(pmf_S2, n_points * n_bins)

        M_ALLOCATE(neighbor_mask_S1_f(n_neighbors, n_points))
        neighbor_mask_S1_f = neighbor_mask_S1
        M_ALLOCATE(neighbor_mask_S2_f(n_neighbors, n_points))
        neighbor_mask_S2_f = neighbor_mask_S2

        call fjct_compute_masked_jsd(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            neighbor_mask_S1 = neighbor_mask_S1_f,&
            neighbor_mask_S2 = neighbor_mask_S2_f,&
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
            ierr = ierr&
        )
    end subroutine fjct_compute_masked_jsd_c

    !> summary: C-wrapper for [[tox_data_integration_per_family(module):fjct_compute_masked_jsd_expert(subroutine)]]
    !| Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts the
    !| residual samples to the neighbors selected by `neighbor_mask_S1`/`neighbor_mask_S2`. Typically
    !| those are all neighbors belonging to one gene family, which is what `fjct_compute_jsd` builds
    !| the masks for from a family index.
    subroutine fjct_compute_masked_jsd_expert_c(&
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
        ) bind(C, name="fjct_compute_masked_jsd_expert_c")
        use tox_data_integration_per_family, only: fjct_compute_masked_jsd_expert

        integer(c_int), intent(in), target :: n_reps_S1
            !! Number of replicates in study 1
        integer(c_int), intent(in), target :: n_reps_S2
            !! Number of replicates in study 2
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors in the studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), target :: neighbor_mask_S1
            !! Mask selecting the neighbors of study 1 to include
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), target :: neighbor_mask_S2
            !! Mask selecting the neighbors of study 2 to include
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_points), intent(out), target :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
        integer(c_int), intent(out), target :: total_included_n_reps
            !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        real(c_double), intent(out), target :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(c_double), dimension(n_points), intent(out), target :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        real(c_double), dimension(n_points, n_bins), intent(out), target :: pmf_S1
            !! Normalized histogram counts for study 1
        real(c_double), dimension(n_points, n_bins), intent(out), target :: pmf_S2
            !! Normalized histogram counts for study 2
        integer(c_int), dimension(n_points, n_bins), intent(out), target :: tmp_counts
            !! Work array for the histogram counts
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(:, :), allocatable :: neighbor_mask_S1_f
        logical, dimension(:, :), allocatable :: neighbor_mask_S2_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_reps_S1)
        M_CHECK_NON_NULL(n_reps_S2)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(total_included_n_reps)
        M_CHECK_NON_NULL(global_js_divergence)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighbor_mask_S1, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighbor_mask_S2, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S1, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S2, n_points)
        M_CHECK_ARRAY_NON_NULL(weights, n_points)
        M_CHECK_ARRAY_NON_NULL(pmf_S1, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(pmf_S2, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, n_points * n_bins)

        M_ALLOCATE(neighbor_mask_S1_f(n_neighbors, n_points))
        neighbor_mask_S1_f = neighbor_mask_S1
        M_ALLOCATE(neighbor_mask_S2_f(n_neighbors, n_points))
        neighbor_mask_S2_f = neighbor_mask_S2

        call fjct_compute_masked_jsd_expert(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            neighbor_mask_S1 = neighbor_mask_S1_f,&
            neighbor_mask_S2 = neighbor_mask_S2_f,&
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
            tmp_counts = tmp_counts,&
            ierr = ierr&
        )
    end subroutine fjct_compute_masked_jsd_expert_c

    !> summary: C-wrapper for [[tox_data_integration_per_family(module):fjct_compute_contribution_scores(subroutine)]]
    !| Combines
    !|
    !| 1. how divergent the family is between the studies (`global_js_divergences`), and
    !| 2. how much residual support the family has overall (`total_included_n_reps_per_f`),
    !|
    !| using the outputs of `fjct_compute_jsd`, collected for the analyzed sub-neighborhoods.
    subroutine fjct_compute_contribution_scores_c(&
            global_js_divergences,&
            total_included_n_reps_per_f,&
            k_families,&
            support_weights,&
            contribution_scores,&
            ierr&
        ) bind(C, name="fjct_compute_contribution_scores_c")
        use tox_data_integration_per_family, only: fjct_compute_contribution_scores

        integer(c_int), intent(in), target :: k_families
            !! Number of sub-neighborhoods analyzed
        real(c_double), dimension(k_families), intent(in), target :: global_js_divergences
            !! Per-sub-neighborhood weighted global JSD
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(k_families), intent(in), target :: total_included_n_reps_per_f
            !! Per-sub-neighborhood `total_included_n_reps`
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(k_families), intent(out), target :: support_weights
            !! Per-sub-neighborhood calculated support weight (ratio between its `total_included_n_reps` and `sum(total_included_n_reps_per_f)`, zero if there were no replicates included at all)
        real(c_double), dimension(k_families), intent(out), target :: contribution_scores
            !! Per-sub-neighborhood calculated contribution ( \( support\_weights_i * global\_js\_divergences_i \) )
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(k_families)
        M_CHECK_ARRAY_NON_NULL(global_js_divergences, k_families)
        M_CHECK_ARRAY_NON_NULL(total_included_n_reps_per_f, k_families)
        M_CHECK_ARRAY_NON_NULL(support_weights, k_families)
        M_CHECK_ARRAY_NON_NULL(contribution_scores, k_families)

        call fjct_compute_contribution_scores(&
            global_js_divergences = global_js_divergences,&
            total_included_n_reps_per_f = total_included_n_reps_per_f,&
            k_families = k_families,&
            support_weights = support_weights,&
            contribution_scores = contribution_scores,&
            ierr = ierr&
        )
    end subroutine fjct_compute_contribution_scores_c

end module tox_data_integration_per_family_c
#endif
