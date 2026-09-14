#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation per family
!|
!| The JSD value for a sub-neighborhood -- typically the genes of one family -- obtained by
!| driving the same pipeline over a masked set of neighbors. Answers whether two studies are
!| compatible *for this family*, which the global figure can hide either way.
module tox_data_integration_per_family_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use f42_math_impl, only: is_close
    use tox_data_integration_jsd_impl, only: jct_compute_jsd_pipeline_helper
    M_IMPLICIT_NONE

contains

    !> summary: Compute the family-level compatibility score between two studies for a single gene family
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts
    !| the residual samples to the genes belonging to the family `family_idx`.
    pure subroutine fjct_compute_jsd_impl(family_idx, gene_to_family_S1, gene_to_family_S2, n_genes_S1, n_genes_S2, neighborhood_residuals_S1, neighborhood_residuals_S2, &
                                            neighborhood_genes_S1, neighborhood_genes_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, js_divergences, &
                                            included_n_reps_S1, included_n_reps_S2, total_included_n_reps, global_js_divergence, weights, &
                                            tmp_neighbor_mask_S1, tmp_neighbor_mask_S2, tmp_pmf_S1, tmp_pmf_S2, tmp_counts &
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
            !! DM_MIN(1_int32)
        integer(int32), dimension(n_genes_S1), intent(in) :: gene_to_family_S1
            !! Mapping for study 1: Each index (gene) holds the index of its family
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_genes_S2), intent(in) :: gene_to_family_S2
            !! Mapping for study 2: Each index (gene) holds the index of its family
            !! DM_MIN(0_int32)
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! DM_ALLOW_NAN
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! DM_ALLOW_NAN
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S1
            !! Indices of the selected neighborhood genes of study 1
            !! DM_MIN(1_int32)
            !! DM_MAX(n_genes_S1)
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S2
            !! Indices of the selected neighborhood genes of study 2
            !! DM_MIN(1_int32)
            !! DM_MAX(n_genes_S2)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies
            !! DM_MIN(0.0_real64)
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

        integer(int32) :: i_point, i_neighbor

        ! Set up mask for filtered analysis -> only include neighbors being part of the family
        do concurrent(i_point=1:n_points)
            do concurrent(i_neighbor=1:n_neighbors) shared(family_idx, tmp_neighbor_mask_S1, tmp_neighbor_mask_S2, gene_to_family_S1, neighborhood_genes_S1, gene_to_family_S2, neighborhood_genes_S2)
                tmp_neighbor_mask_S1(i_neighbor, i_point) = gene_to_family_S1(neighborhood_genes_S1(i_neighbor, i_point)) == family_idx
                tmp_neighbor_mask_S2(i_neighbor, i_point) = gene_to_family_S2(neighborhood_genes_S2(i_neighbor, i_point)) == family_idx
            end do
        end do

        call fjct_compute_masked_jsd_impl(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, tmp_neighbor_mask_S1, tmp_neighbor_mask_S2, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, total_included_n_reps, global_js_divergence, weights, tmp_pmf_S1, tmp_pmf_S2, tmp_counts)
    end subroutine fjct_compute_jsd_impl

    !> summary: Compute the compatibility score between two studies for a single masked sub-neighborhood
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts the
    !| residual samples to the neighbors selected by `neighbor_mask_S1`/`neighbor_mask_S2`. Typically
    !| those are all neighbors belonging to one gene family, which is what `fjct_compute_jsd` builds
    !| the masks for from a family index.
    pure subroutine fjct_compute_masked_jsd_impl(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, neighbor_mask_S1, neighbor_mask_S2, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, total_included_n_reps, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts)
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! DM_ALLOW_NAN
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! DM_ALLOW_NAN
        logical(c_bool), dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S1
            !! Mask selecting the neighbors of study 1 to include
        logical(c_bool), dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S2
            !! Mask selecting the neighbors of study 2 to include
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies
            !! DM_MIN(0.0_real64)
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

        call jct_compute_jsd_pipeline_helper(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, neighbor_mask_S1, neighbor_mask_S2)
        total_included_n_reps = sum(included_n_reps_S1) + sum(included_n_reps_S2)
    end subroutine fjct_compute_masked_jsd_impl

    !> summary: Compute the per-sub-neighborhood contribution score
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Combines
    !|
    !| 1. how divergent the family is between the studies (`global_js_divergences`), and
    !| 2. how much residual support the family has overall (`total_included_n_reps_per_f`),
    !|
    !| using the outputs of `fjct_compute_jsd`, collected for the analyzed sub-neighborhoods.
    pure subroutine fjct_compute_contribution_scores_impl(global_js_divergences, total_included_n_reps_per_f, k_families, support_weights, contribution_scores)
        integer(int32), intent(in) :: k_families
            !! Number of sub-neighborhoods analyzed
        integer(int32), dimension(k_families), intent(in) :: total_included_n_reps_per_f
            !! Per-sub-neighborhood `total_included_n_reps`
            !! DM_MIN(0_int32)
        real(real64), dimension(k_families), intent(in) :: global_js_divergences
            !! Per-sub-neighborhood weighted global JSD
            !! DM_MIN(0.0_real64)
        real(real64), dimension(k_families), intent(out) :: support_weights
            !! Per-sub-neighborhood calculated support weight (ratio between its `total_included_n_reps` and `sum(total_included_n_reps_per_f)`, zero if there were no replicates included at all)
        real(real64), dimension(k_families), intent(out) :: contribution_scores
            !! Per-sub-neighborhood calculated contribution ( \( support\_weights_i * global\_js\_divergences_i \) )

        integer(int32) :: i_family
        real(real64) :: total_included_n_reps

        total_included_n_reps = real(sum(total_included_n_reps_per_f), kind=real64)

        ! Guard against a zero (or near-zero) total support: without this early `return`, execution would
        ! fall through into the loop below and divide by `total_included_n_reps`, producing Inf/NaN scores
        ! for a family with no included replicates instead of the well-defined zero score.
        if (is_close(total_included_n_reps, 0.0_real64)) then
            support_weights = 0.0_real64
            contribution_scores = 0.0_real64
            return
        end if

        do concurrent(i_family=1:k_families) shared(support_weights, total_included_n_reps_per_f, total_included_n_reps, contribution_scores, global_js_divergences)
            support_weights(i_family) = real(total_included_n_reps_per_f(i_family), kind=real64)/total_included_n_reps

            contribution_scores(i_family) = support_weights(i_family)*global_js_divergences(i_family)
        end do
    end subroutine fjct_compute_contribution_scores_impl

end module tox_data_integration_per_family_impl
