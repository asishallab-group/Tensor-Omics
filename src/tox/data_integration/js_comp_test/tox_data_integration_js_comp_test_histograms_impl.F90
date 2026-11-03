#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Histograms
!|
!| The per-study residual histograms of the JS-Comp-Test, built at the bin counts and residual
!| ranges the bin-count search
!| ([[tox_data_integration_js_comp_test_binning_impl(module)]]) chose per reference point, by
!| [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]: on fixed-k
!| neighborhoods every study at once, on adaptive neighborhoods one (reference point, study) pair at
!| a time. Nothing here is published on its own; the final runs and the parameter searches of
!| [[tox_data_integration_js_comp_test_run_impl(module)]],
!| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module)]] and
!| [[tox_data_integration_js_comp_test_adaptive_search_impl(module)]] build their histograms with it.
module tox_data_integration_js_comp_test_histograms_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_data_integration_jsd_impl, only: build_residual_histograms_impl
    use tox_data_integration_preprocessing_impl, only: materialize_pooled_neighborhood
    use tox_data_integration_js_comp_test_binning_impl, only: MAX_N_BINS
    M_IMPLICIT_NONE
    private
    public :: build_point_study_histogram, build_fixed_k_histograms, build_adaptive_histograms

contains

    !> summary: Build one study's residual histogram for one reference point
    !| AUTHOR_LASZLO_LANG
    !| For ONE (reference point, study) pair: gathers the neighbors' residuals and bins them with
    !| [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]] for a
    !| single point, which gives exactly the matching row of a call over all points. Only row
    !| `i_point` of `counts`/`pmf` is written, zero beyond `n_bins`; NaN residuals are skipped.
    !| Internal helper, not published.
    !|
    !| Precondition, not checked here: every entry of `neighbor_indices` is in
    !| `[1, max_n_genes_all_studies]`, `1 <= n_bins <= max_n_bins <= MAX_N_BINS`, and
    !| `1 <= i_point <= n_points`.
    pure subroutine build_point_study_histogram(residuals_study, max_n_reps_all_studies, max_n_genes_all_studies, &
                                                n_neighbors, neighbor_indices, shared_residual_range_low, &
                                                shared_residual_range_high, n_bins, n_points, max_n_bins, i_point, &
                                                counts, pmf, included_n_reps, tmp_neighbor_residuals)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors of this reference point in this study
        real(real64), intent(in) :: residuals_study(max_n_reps_all_studies, max_n_genes_all_studies)
            !! This study's signed residuals, NaN allowed for missing values
        integer(int32), intent(in) :: neighbor_indices(n_neighbors)
            !! Gene indices of this reference point's neighbors in this study
        real(real64), intent(in) :: shared_residual_range_low
            !! Lower bound of this reference point's histogram range (R_low)
        real(real64), intent(in) :: shared_residual_range_high
            !! Upper bound of this reference point's histogram range (R_high)
        integer(int32), intent(in) :: n_bins
            !! This reference point's bin count
        integer(int32), intent(in) :: n_points
            !! First extent of `counts`/`pmf`: the number of reference points
        integer(int32), intent(in) :: max_n_bins
            !! Second extent of `counts`/`pmf`: the widest bin count of any reference point
        integer(int32), intent(in) :: i_point
            !! The row of `counts`/`pmf` this reference point's histogram is written to
        integer(int32), intent(inout) :: counts(n_points, max_n_bins)
            !! The study's histogram counts; row `i_point` is overwritten, zero beyond `n_bins`
        real(real64), intent(inout) :: pmf(n_points, max_n_bins)
            !! The study's pmfs; row `i_point` is overwritten (all zero when `included_n_reps` is 0)
        integer(int32), intent(out) :: included_n_reps
            !! Number of non-NaN residuals that were binned
        real(real64), intent(out) :: tmp_neighbor_residuals(max_n_reps_all_studies, n_neighbors)
            !! Working array: the gathered residuals of this reference point's neighbors

        integer(int32) :: i_neighbor
        integer(int32) :: point_n_bins(1), point_included_n_reps(1)
        integer(int32) :: point_counts(1, MAX_N_BINS)
        real(real64) :: point_range_low(1), point_range_high(1)
        real(real64) :: point_pmf(1, MAX_N_BINS)

        do concurrent(i_neighbor=1:n_neighbors) shared(tmp_neighbor_residuals, residuals_study, neighbor_indices)
            tmp_neighbor_residuals(:, i_neighbor) = residuals_study(:, neighbor_indices(i_neighbor))
        end do

        point_n_bins(1) = n_bins
        point_range_low(1) = shared_residual_range_low
        point_range_high(1) = shared_residual_range_high

        ! point_counts/point_pmf passed whole: with one row, the callee's (1, max_n_bins) dummies are
        ! exactly their first max_n_bins elements (max_n_bins <= MAX_N_BINS by the precondition).
        call build_residual_histograms_impl(tmp_neighbor_residuals, max_n_reps_all_studies, n_neighbors, 1_int32, &
                                            point_range_low, point_range_high, max_n_bins, point_n_bins, &
                                            point_counts, point_pmf, &
                                            point_included_n_reps)

        counts(i_point, :) = point_counts(1, 1:max_n_bins)
        pmf(i_point, :) = point_pmf(1, 1:max_n_bins)
        included_n_reps = point_included_n_reps(1)
    end subroutine build_point_study_histogram

    !> summary: Build one study's fixed-k histograms for every reference point
    !| AUTHOR_LASZLO_LANG
    !| Pass C of the fixed-k pipeline for ONE study: gathers every reference point's neighbor
    !| residuals into `tmp_neighborhood_residuals_gathered`, then builds all of the study's
    !| histograms with one
    !| [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]] call at
    !| the per-point bin counts. Private helper of `build_fixed_k_histograms`. Precondition, not
    !| checked here: every index is in `[1, max_n_genes_all_studies]` and
    !| `1 <= n_bins_per_point <= max_n_bins`.
    pure subroutine build_fixed_k_study_histograms(residuals_study, max_n_reps_all_studies, max_n_genes_all_studies, &
                                                   n_neighbors, n_points, neighborhood_indices_study, &
                                                   shared_residual_range_low, shared_residual_range_high, &
                                                   n_bins_per_point, max_n_bins, counts, pmfs, included_n_reps, &
                                                   tmp_neighborhood_residuals_gathered)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per reference point
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        real(real64), intent(in) :: residuals_study(max_n_reps_all_studies, max_n_genes_all_studies)
            !! This study's signed residuals, NaN allowed for missing values
        integer(int32), intent(in) :: neighborhood_indices_study(n_neighbors, n_points)
            !! This study's neighbor gene indices per reference point
        real(real64), intent(in) :: shared_residual_range_low(n_points)
            !! Lower histogram range bound (R_low) per reference point
        real(real64), intent(in) :: shared_residual_range_high(n_points)
            !! Upper histogram range bound (R_high) per reference point
        integer(int32), intent(in) :: n_bins_per_point(n_points)
            !! Bin count per reference point
        integer(int32), intent(in) :: max_n_bins
            !! Second extent of `counts`/`pmfs`: the widest bin count of any reference point
        integer(int32), intent(out) :: counts(n_points, max_n_bins)
            !! Residual count per bin, per reference point; zero beyond each point's bin count
        real(real64), intent(out) :: pmfs(n_points, max_n_bins)
            !! `counts` normalized per reference point
        integer(int32), intent(out) :: included_n_reps(n_points)
            !! Number of non-NaN residuals binned per reference point
        real(real64), intent(out) :: tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, n_neighbors, n_points)
            !! Working array: the gathered residuals of every reference point's neighbors

        integer(int32) :: i_point, i_neighbor, gene_idx

        do concurrent(i_point=1:n_points, i_neighbor=1:n_neighbors) local(gene_idx) &
                shared(tmp_neighborhood_residuals_gathered, residuals_study, neighborhood_indices_study)
            gene_idx = neighborhood_indices_study(i_neighbor, i_point)
            tmp_neighborhood_residuals_gathered(:, i_neighbor, i_point) = residuals_study(:, gene_idx)
        end do

        call build_residual_histograms_impl(tmp_neighborhood_residuals_gathered, max_n_reps_all_studies, n_neighbors, &
                                            n_points, shared_residual_range_low, shared_residual_range_high, max_n_bins, &
                                            n_bins_per_point, counts, pmfs, included_n_reps)
    end subroutine build_fixed_k_study_histograms

    !> summary: Build every (reference point, study) histogram on adaptive neighborhoods
    !| AUTHOR_LASZLO_LANG
    !| Pass C on adaptive neighborhoods: decodes every reference point's pooled range into
    !| per-study gene lists
    !| ([[tox_data_integration_preprocessing_impl(module):materialize_pooled_neighborhood(subroutine)]])
    !| and builds each (point, study) histogram with `build_point_study_histogram` into bins
    !| `1:max_n_bins` of `counts`/`pmfs`. Sequential, so one gather buffer serves every pair.
    !| Internal helper, not published; public because `run_js_comp_test_adaptive_impl`
    !| (`n_bins_extent = M_MAX_N_BINS`) and the adaptive parameter search (`n_bins_extent =
    !| max_n_bins`, its buffers viewed packed) call it.
    pure subroutine build_adaptive_histograms(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                              n_points, max_n_bins, n_bins_extent, gene_means_perm_all, &
                                              pooled_neighborhood_range, n_neighbors_per_point, &
                                              shared_residual_range_low, shared_residual_range_high, n_bins_per_point, &
                                              counts, pmfs, included_n_reps, tmp_point_neighborhood_indices, &
                                              tmp_neighbor_residuals)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: max_n_bins
            !! Widest per-point bin count; only bins `1:max_n_bins` of `counts`/`pmfs` are written
        integer(int32), intent(in) :: n_bins_extent
            !! Second extent of `counts`/`pmfs`, at least `max_n_bins`
        real(real64), intent(in) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Matrix of signed residuals per study, NaN allowed for missing values
        integer(int32), intent(in) :: gene_means_perm_all(max_n_genes_all_studies*n_studies)
            !! Sorting permutation of the pooled gene means, as the construction sorts them
        integer(int32), intent(in) :: pooled_neighborhood_range(2, n_points)
            !! Each reference point's `[first, last]` pooled position
        integer(int32), intent(out) :: n_neighbors_per_point(n_studies, n_points)
            !! Each reference point's per-study gene count, as decoded
        real(real64), intent(in) :: shared_residual_range_low(n_points)
            !! Lower histogram range bound (R_low) per reference point
        real(real64), intent(in) :: shared_residual_range_high(n_points)
            !! Upper histogram range bound (R_high) per reference point
        integer(int32), intent(in) :: n_bins_per_point(n_points)
            !! Bin count per reference point
        integer(int32), intent(inout) :: counts(n_points, n_bins_extent, n_studies)
            !! Every study's histogram counts; bins `1:max_n_bins` are overwritten
        real(real64), intent(inout) :: pmfs(n_points, n_bins_extent, n_studies)
            !! Every study's pmfs; bins `1:max_n_bins` are overwritten
        integer(int32), intent(out) :: included_n_reps(n_points, n_studies)
            !! Number of non-NaN residuals binned per reference point, per study
        integer(int32), intent(out) :: tmp_point_neighborhood_indices(max_n_genes_all_studies, n_studies)
            !! Working array: one reference point's gene indices per study, reused per point
        real(real64), intent(out) :: tmp_neighbor_residuals(max_n_reps_all_studies, max_n_genes_all_studies)
            !! Working array: one (reference point, study) pair's gathered residuals, reused per pair

        integer(int32) :: i_point, i_study, n_study_neighbors

        do i_point = 1, n_points
            call materialize_pooled_neighborhood(gene_means_perm_all, max_n_genes_all_studies, n_studies, &
                                                 pooled_neighborhood_range(1, i_point), &
                                                 pooled_neighborhood_range(2, i_point), max_n_genes_all_studies, &
                                                 tmp_point_neighborhood_indices, n_neighbors_per_point(:, i_point))
            do i_study = 1, n_studies
                n_study_neighbors = n_neighbors_per_point(i_study, i_point)
                call build_point_study_histogram(residuals(:, :, i_study), max_n_reps_all_studies, &
                                                 max_n_genes_all_studies, n_study_neighbors, &
                                                 tmp_point_neighborhood_indices(1:n_study_neighbors, i_study), &
                                                 shared_residual_range_low(i_point), shared_residual_range_high(i_point), &
                                                 n_bins_per_point(i_point), n_points, max_n_bins, i_point, &
                                                 counts(:, 1:max_n_bins, i_study), pmfs(:, 1:max_n_bins, i_study), &
                                                 included_n_reps(i_point, i_study), &
                                                 tmp_neighbor_residuals(:, 1:n_study_neighbors))
            end do
        end do
    end subroutine build_adaptive_histograms

    !> summary: Build every study's fixed-k histograms for every reference point
    !| AUTHOR_LASZLO_LANG
    !| Pass C of the fixed-k pipeline: every study's histograms, one
    !| `build_fixed_k_study_histograms` call per study, into bins `1:max_n_bins` of `counts`/`pmfs`.
    !| Sequential, so one gather buffer serves every study. Internal helper, not published; public
    !| because `run_js_comp_test_impl` (`n_bins_extent = M_MAX_N_BINS`) and the fixed-k parameter
    !| search (`n_bins_extent = max_n_bins`, its buffers viewed packed at the candidate's sizes) call
    !| it. Precondition, not checked here: every index is in `[1, max_n_genes_all_studies]` and
    !| `1 <= n_bins_per_point <= max_n_bins <= n_bins_extent`.
    pure subroutine build_fixed_k_histograms(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                             n_neighbors, n_points, max_n_bins, n_bins_extent, neighborhood_indices, &
                                             shared_residual_range_low, shared_residual_range_high, n_bins_per_point, &
                                             counts, pmfs, included_n_reps, tmp_neighborhood_residuals_gathered)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per reference point
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: max_n_bins
            !! Widest per-point bin count; only bins `1:max_n_bins` of `counts`/`pmfs` are written
        integer(int32), intent(in) :: n_bins_extent
            !! Second extent of `counts`/`pmfs`, at least `max_n_bins`
        real(real64), intent(in) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Matrix of signed residuals per study, NaN allowed for missing values
        integer(int32), intent(in) :: neighborhood_indices(n_neighbors, n_points, n_studies)
            !! Every study's neighbor gene indices per reference point
        real(real64), intent(in) :: shared_residual_range_low(n_points)
            !! Lower histogram range bound (R_low) per reference point
        real(real64), intent(in) :: shared_residual_range_high(n_points)
            !! Upper histogram range bound (R_high) per reference point
        integer(int32), intent(in) :: n_bins_per_point(n_points)
            !! Bin count per reference point
        integer(int32), intent(inout) :: counts(n_points, n_bins_extent, n_studies)
            !! Every study's histogram counts; bins `1:max_n_bins` are overwritten
        real(real64), intent(inout) :: pmfs(n_points, n_bins_extent, n_studies)
            !! Every study's pmfs; bins `1:max_n_bins` are overwritten
        integer(int32), intent(out) :: included_n_reps(n_points, n_studies)
            !! Number of non-NaN residuals binned per reference point, per study
        real(real64), intent(out) :: tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, n_neighbors, n_points)
            !! Working array: one study's gathered residuals, reused per study

        integer(int32) :: i_study

        do i_study = 1, n_studies
            call build_fixed_k_study_histograms(residuals(:, :, i_study), max_n_reps_all_studies, max_n_genes_all_studies, &
                                                n_neighbors, n_points, neighborhood_indices(:, :, i_study), &
                                                shared_residual_range_low, shared_residual_range_high, n_bins_per_point, &
                                                max_n_bins, counts(:, 1:max_n_bins, i_study), pmfs(:, 1:max_n_bins, i_study), &
                                                included_n_reps(:, i_study), tmp_neighborhood_residuals_gathered)
        end do
    end subroutine build_fixed_k_histograms

end module tox_data_integration_js_comp_test_histograms_impl
