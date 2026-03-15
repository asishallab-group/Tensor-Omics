#include "macros.h"

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) - JSD Calculation
!|
!| This module implements the pipeline to obtain the JSD value from neighborhood residuals obtained from [[tox_data_integration_preprocessing(module)]].
module tox_data_integration_jsd
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_positive_inf
    use, intrinsic :: iso_c_binding, only: c_ptr
    use f42_utils, only: clamp, is_close, LOG_2
    use tox_errors, only: set_ok, set_err, is_err, validate_dimension_size, validate_in_range_real, validate_all_in_range_real, validate_all_in_range_int
    implicit none

contains

    !> Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
    pure subroutine build_residual_histograms(neighborhood_residuals, n_neighbors, n_points, residuals, max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins, counts, pmf, included_n_reps, ierr, neighbor_mask)
        integer(int32), intent(in) :: n_neighbors
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_residuals
            !! Computed neighborhood residuals for a study ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(in) :: residuals
            !! Matrix of signed residuals for one study
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R) from [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        integer(int32), dimension(n_bins, n_points), intent(out) :: counts
            !! Absolute counts of a residual per bin
        real(real64), dimension(n_bins, n_points), intent(out) :: pmf
            !! `counts` normalized to `0 <= counts(i, :) <= 1` and `sum(counts(i, :)) == 1`
        integer(int32), dimension(n_points), intent(out) :: included_n_reps
            !! Stores the count of non-NaN replicates (included ones)
        integer(int32), intent(out) :: ierr
            !! Error code
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
            !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)

        call set_ok(ierr)

        call validate_dimension_size(n_neighbors, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_dimension_size(max_n_reps_all_studies, ierr, arg_pos=5_int32)
        call validate_dimension_size(max_n_genes_all_studies, ierr, arg_pos=6_int32)
        call validate_in_range_real(shared_residual_range, ierr, min=0.0_real64, arg_pos=7_int32)
        call validate_dimension_size(n_bins, ierr, arg_pos=8_int32)
        call validate_all_in_range_int(neighborhood_residuals, size(neighborhood_residuals, kind=int32), ierr, min=1_int32, max=max_n_genes_all_studies, arg_pos=1_int32)

        if (is_err(ierr)) return

        call build_residual_histograms_helper(neighborhood_residuals, n_neighbors, n_points, residuals, max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins, counts, pmf, included_n_reps, neighbor_mask)
    end subroutine build_residual_histograms

    !> (no input validation) Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
    pure subroutine build_residual_histograms_helper(neighborhood_residuals, n_neighbors, n_points, residuals, max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins, counts, pmf, included_n_reps, neighbor_mask)
        integer(int32), intent(in) :: n_neighbors
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_residuals
            !! Computed neighborhood residuals for a study ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(in) :: residuals
            !! Matrix of signed residuals for one study
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R) from [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        integer(int32), dimension(n_bins, n_points), intent(out) :: counts
            !! Absolute counts of a residual per bin
        real(real64), dimension(n_bins, n_points), intent(out) :: pmf
            !! `counts` normalized to `0 <= counts(i, :) <= 1` and `sum(counts(i, :)) == 1`
        integer(int32), dimension(n_points), intent(out) :: included_n_reps
            !! Stores the count of non-NaN replicates (included ones)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
            !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)

        real(real64) :: bin_width, clamped_residual, replicate
        integer(int32) :: bin_idx, i_neighbor, i_rep, included_reps, i_point
        logical :: filter_neighbors

        bin_width = 2.0_real64 * shared_residual_range / real(n_bins, real64)
        counts = 0_int32
        pmf = 0.0_real64

        filter_neighbors = present(neighbor_mask)

        ! 1. assign the bins to the residuals (increase the respective count)
        ! outer loop cannot be concurrent, as counts of same residuals and bins but different neighbors might be changed at the same time
        do concurrent (i_point = 1:n_points) local(included_reps) shared(included_n_reps)
            included_reps = 0_int32
            do concurrent (i_neighbor = 1:n_neighbors) &
                    local(i_rep, clamped_residual, bin_idx, replicate) &
                    shared(filter_neighbors, max_n_reps_all_studies, counts, neighborhood_residuals, shared_residual_range, bin_width) &
                    reduce(+:included_reps)
                ! Exclude neighbor if desired
                if (filter_neighbors) then
                    if (.not. neighbor_mask(i_neighbor, i_point)) cycle
                end if

                ! Count non-NaNs and assign the to a bin
                do i_rep = 1, max_n_reps_all_studies
                    replicate = residuals(i_rep, neighborhood_residuals(i_neighbor, i_point))
                    if (.not. ieee_is_nan(replicate)) then
                        ! clamp residual to histogram range
                        clamped_residual = clamp(replicate, min_val=-shared_residual_range, max_val=shared_residual_range)

                        ! assign bin to residual
                        bin_idx = min(n_bins, int( (clamped_residual + shared_residual_range) / bin_width ) + 1)
                        counts(bin_idx, i_point) = counts(bin_idx, i_point) + 1

                        included_reps = included_reps + 1
                    end if
                end do
            end do
            included_n_reps(i_point) = included_reps
        end do

        call calc_pmf_helper(counts, pmf, included_n_reps, n_bins, n_points)
    end subroutine build_residual_histograms_helper

    !> Helper for calculating the pmf for histogram counts (e.g. in [[tox_data_integration(module):build_residual_histograms_helper(interface)]])
    pure subroutine calc_pmf_helper(counts, pmf, included_n_reps, n_bins, n_points)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), dimension(n_bins, n_points), intent(out) :: pmf
            !! `counts` normalized to `0 <= counts(i, :) <= 1` and `sum(counts(i, :)) == 1`
        integer(int32), dimension(n_bins, n_points), intent(in) :: counts
            !! Absolute counts of a residual per bin
        integer(int32), dimension(n_points), intent(in) :: included_n_reps
            !! Stores the count of non-NaN replicates (included ones in `counts`)

        integer(int32) :: i_point, i_bin

        do concurrent (i_point = 1:n_points)
            do concurrent (i_bin = 1:n_bins) shared(pmf, i_point, included_n_reps, counts)
                if (included_n_reps(i_point) == 0) then
                    pmf(i_bin, i_point) = 0.0_real64
                else
                    pmf(i_bin, i_point) = real(counts(i_bin, i_point), real64) / real(included_n_reps(i_point), real64)
                end if
            end do
        end do
    end subroutine calc_pmf_helper

    !> Having the probabilities `pmf` from [[tox_data_integration(module):build_residual_histograms(interface)]], this subroutine computes the Jensen-Shannon divergence per reference point/neighbor
    pure subroutine compute_divergence_per_reference_point(pmf_S1, pmf_S2, n_points, n_bins, js_divergences, ierr)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(real64), dimension(n_bins, n_points), intent(in) :: pmf_S1
            !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 1
        real(real64), dimension(n_bins, n_points), intent(in) :: pmf_S2
            !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 2
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_bins, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(pmf_S1, size(pmf_S1, kind=int32), ierr, arg_pos=1_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(pmf_S2, size(pmf_S2, kind=int32), ierr, arg_pos=2_int32, min=0.0_real64, max=1.0_real64)

        if (is_err(ierr)) return

        call compute_divergence_per_reference_point_helper(pmf_S1, pmf_S2, n_points, n_bins, js_divergences)
    end subroutine compute_divergence_per_reference_point

    !> (no input validation) Having the probabilities `pmf` from [[tox_data_integration(module):build_residual_histograms(interface)]], this subroutine computes the Jensen-Shannon divergence per reference point/neighbor
    pure subroutine compute_divergence_per_reference_point_helper(pmf_S1, pmf_S2, n_points, n_bins, js_divergences)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(real64), dimension(n_bins, n_points), intent(in) :: pmf_S1
            !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 1
        real(real64), dimension(n_bins, n_points), intent(in) :: pmf_S2
            !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 2
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point

        real(real64) :: S_mean, s1_val, s2_val
        integer(int32) :: i_bin, i_point

        js_divergences = 0.0_real64

        ! 1. compute the Knullback-Leibler (KL) divergences
        ! Note that the J-S divergence is defined as `0.5 * KL_S1 + 0.5 * KL_S2`, equivalent to `0.5 * (KL_S1 + KL_S2)`.
        ! Thus, instead of computing KL_S* independently, it accumulates directly in the `js_divergences` output.
        ! Another thing, switching the loops would enable both to run concurrently and the 0.5*js_divergences step could be done in one go,
        ! but cache locality still beats that, except for the case of thousands of neighbors, which might not be the common case.
        do i_bin = 1, n_bins
            do concurrent (i_point = 1:n_points) local(s1_val, s2_val, S_mean) shared(i_bin, pmf_S1, pmf_S2, js_divergences)
                s1_val = pmf_S1(i_bin, i_point)
                s2_val = pmf_S2(i_bin, i_point)
                S_mean = 0.5_real64 * (s1_val + s2_val)

                if (.not. is_close(S_mean, 0.0_real64)) then
                    if (s1_val > 0.0_real64) then
                        js_divergences(i_point) = js_divergences(i_point) + s1_val * log(s1_val / S_mean)
                    end if

                    if (s2_val > 0.0_real64) then
                        js_divergences(i_point) = js_divergences(i_point) + s2_val * log(s2_val / S_mean)
                    end if
                end if
            end do
        end do

        ! 2. Compute the js_divergences
        do concurrent (i_point = 1:n_points) shared(js_divergences)
            js_divergences(i_point) = (0.5_real64 * js_divergences(i_point)) / LOG_2  ! Div by ln2 to have a 0-1 scaling instead of 0-ln2
        end do
    end subroutine compute_divergence_per_reference_point_helper

    !> Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
    pure subroutine compute_weighted_global_divergence(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, ierr)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        real(real64), dimension(n_points), intent(in) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(js_divergences, size(js_divergences, kind=int32), ierr, arg_pos=1_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_int(included_n_reps_S1, size(included_n_reps_S1, kind=int32), ierr, arg_pos=3_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps_S2, size(included_n_reps_S2, kind=int32), ierr, arg_pos=4_int32, min=0_int32)

        if (is_err(ierr)) return

        call compute_weighted_global_divergence_helper(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights)
    end subroutine compute_weighted_global_divergence

    !> (no input validation) Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
    pure subroutine compute_weighted_global_divergence_helper(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        real(real64), dimension(n_points), intent(in) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`

        integer(int32) :: i_point, included_reps
        real(real64) :: total_sample_count

        global_js_divergence = 0.0_real64

        total_sample_count = real(sum(included_n_reps_S1) + sum(included_n_reps_S2), real64)

        if (is_close(total_sample_count, 0.0_real64)) then
            weights = 0.0_real64
        else
            ! Calculate the global Jensen-Shannon divergence
            do concurrent (i_point = 1:n_points) local(included_reps) shared(weights, included_n_reps_S1, included_n_reps_S2, total_sample_count) reduce(+:global_js_divergence)
                included_reps = included_n_reps_S1(i_point) + included_n_reps_S2(i_point)

                weights(i_point) = real(included_reps, real64) / total_sample_count

                global_js_divergence = global_js_divergence + weights(i_point) * js_divergences(i_point)
            end do
        end if
    end subroutine compute_weighted_global_divergence_helper
end module tox_data_integration_jsd