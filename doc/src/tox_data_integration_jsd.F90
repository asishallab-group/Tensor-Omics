#include "macros.h"

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation
!|
!| This module implements the pipeline to obtain the JSD value from neighborhood residuals obtained from [[tox_data_integration_preprocessing(submodule)]].
submodule (tox_data_integration) tox_data_integration_jsd
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use f42_utils, only: clamp, calc_percentile_helper, is_close, sort_array_heapsort, shuffle_vector, init_random
    use tox_errors, only: set_ok, set_err, is_err, ERR_ALLOC_FAIL, validate_dimension_size, validate_in_range_real, validate_all_in_range_real, validate_in_range_int, validate_all_in_range_int
    implicit none
contains

    !> Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
    pure module subroutine determine_shared_residual_range(abs_residual_pool, abs_residual_pool_perm, pool_size, shared_residual_range, ierr, residual_range_quantile)
        integer(int32), intent(in) :: pool_size
            !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), dimension(pool_size), intent(in) :: abs_residual_pool
            !! The absolute residual values of the concatenated S1,S2 residuals
        integer(int32), dimension(pool_size), intent(in) :: abs_residual_pool_perm
            !! The permutation vector that sorts `abs_residual_pool`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(pool_size, ierr)
        call validate_in_range_real(residual_range_quantile, ierr, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_int(abs_residual_pool_perm, pool_size, ierr, min=1_int32, max=pool_size)

        if (is_err(ierr)) return

        call determine_shared_residual_range_helper(abs_residual_pool, abs_residual_pool_perm, pool_size, shared_residual_range, residual_range_quantile)
    end subroutine determine_shared_residual_range

    !> (no input validation) Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
    pure module subroutine determine_shared_residual_range_helper(abs_residual_pool, abs_residual_pool_perm, pool_size, shared_residual_range, residual_range_quantile)
        integer(int32), intent(in) :: pool_size
            !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), dimension(pool_size), intent(in) :: abs_residual_pool
            !! The absolute residual values of the concatenated S1,S2 residuals
        integer(int32), dimension(pool_size), intent(in) :: abs_residual_pool_perm
            !! The permutation vector that sorts `abs_residual_pool`

        integer(int32) :: i_pool, last_non_nan
        real(real64) :: actual_quantile

        M_DEFAULT_VAL(residual_range_quantile, actual_quantile, 95.0_real64)

        last_non_nan = pool_size
        ! NaN is always last -> find last non-NaN index for percentile calculation
        do i_pool = last_non_nan, 1, -1
            if (ieee_is_nan(abs_residual_pool(abs_residual_pool_perm(i_pool)))) then
                last_non_nan = last_non_nan - 1
            else
                exit
            end if
        end do

        if (last_non_nan == 0) then
            shared_residual_range = 0.0_real64
            return
        end if

        call calc_percentile_helper(abs_residual_pool(:last_non_nan), abs_residual_pool_perm(:last_non_nan), actual_quantile, shared_residual_range)
    end subroutine determine_shared_residual_range_helper

    !> Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
    pure module subroutine determine_shared_residual_range_alloc(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, shared_residual_range, ierr, residual_range_quantile)
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
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_rep, i_neighbor, pool_size, pool_idx, n_predecessors, i_point
        real(real64), dimension(:, :, :), allocatable, target :: abs_residual_pool
        real(real64), dimension(:), pointer :: abs_residual_pool_flat
        integer(int32), dimension(:, :, :), allocatable, target :: perm
        integer(int32), dimension(:), pointer :: perm_flat

        call set_ok(ierr)

        call validate_dimension_size(n_reps_S1, ierr)
        call validate_dimension_size(n_reps_S2, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_points, ierr)

        if (is_err(ierr)) return

        M_ALLOCATE(abs_residual_pool(n_reps_S1 + n_reps_S2, n_neighbors, n_points))
        M_ALLOCATE(perm(n_reps_S1 + n_reps_S2, n_neighbors, n_points))

        pool_size = size(abs_residual_pool, kind=int32)

        ! Collect the absolute residual values and compute the percentile value
        do concurrent (i_point = 1:n_points)
            do concurrent (i_neighbor = 1:n_neighbors) local(n_predecessors) shared(i_point, n_reps_S1, n_reps_S2)
                n_predecessors = ((i_point - 1) * n_neighbors + (i_neighbor - 1)) * (n_reps_S1 + n_reps_S2)
                do concurrent (i_rep = 1:n_reps_S1) shared(n_predecessors, perm, neighborhood_residuals_S1, abs_residual_pool, i_point, i_neighbor)
                    abs_residual_pool(i_rep, i_neighbor, i_point) = abs(neighborhood_residuals_S1(i_rep, i_neighbor, i_point))

                    perm(i_rep, i_neighbor, i_point) = n_predecessors + i_rep
                end do

                do concurrent (i_rep = 1:n_reps_S2) local(pool_idx) shared(n_predecessors, perm, neighborhood_residuals_S2, abs_residual_pool, i_point, i_neighbor)
                    pool_idx = i_rep + n_reps_S1

                    abs_residual_pool(pool_idx, i_neighbor, i_point) = abs(neighborhood_residuals_S2(i_rep, i_neighbor, i_point))

                    perm(pool_idx, i_neighbor, i_point) = n_predecessors + pool_idx
                end do
            end do
        end do

        abs_residual_pool_flat(1:pool_size) => abs_residual_pool
        perm_flat(1:pool_size) => perm

        call sort_array_heapsort(abs_residual_pool_flat, perm_flat)

        call determine_shared_residual_range(abs_residual_pool_flat, perm_flat, pool_size, shared_residual_range, ierr, residual_range_quantile)
    end subroutine determine_shared_residual_range_alloc

    !> Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
    pure module subroutine build_residual_histograms(neighborhood_residuals, n_reps, n_neighbors, n_points, shared_residual_range, n_bins, counts, pmf, included_n_reps, ierr, neighbor_mask)
        integer(int32), intent(in) :: n_reps
            !! Number of replicates of the study
        integer(int32), intent(in) :: n_neighbors
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps, n_neighbors, n_points), intent(in) :: neighborhood_residuals
            !! Computed neighborhood residuals for a study ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R) from [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        integer(int32), dimension(n_points, n_bins), intent(out) :: counts
            !! Absolute counts of a residual per bin
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf
            !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
        integer(int32), dimension(n_points), intent(out) :: included_n_reps
            !! Stores the count of non-NaN replicates (included ones)
        integer(int32), intent(out) :: ierr
            !! Error code
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
            !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)

        call set_ok(ierr)

        call validate_dimension_size(n_reps, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_dimension_size(n_bins, ierr)
        call validate_in_range_real(shared_residual_range, ierr, min=0.0_real64)

        if (is_err(ierr)) return

        call build_residual_histograms_helper(neighborhood_residuals, n_reps, n_neighbors, n_points, shared_residual_range, n_bins, counts, pmf, included_n_reps, neighbor_mask)
    end subroutine build_residual_histograms

    !> (no input validation) Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
    pure module subroutine build_residual_histograms_helper(neighborhood_residuals, n_reps, n_neighbors, n_points, shared_residual_range, n_bins, counts, pmf, included_n_reps, neighbor_mask)
        integer(int32), intent(in) :: n_reps
            !! Number of replicates of the study
        integer(int32), intent(in) :: n_neighbors
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the study
        real(real64), dimension(n_reps, n_neighbors, n_points), intent(in) :: neighborhood_residuals
            !! Computed neighborhood residuals for a study ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R) from [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        integer(int32), dimension(n_points, n_bins), intent(out) :: counts
            !! Absolute counts of a residual per bin
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf
            !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
        integer(int32), dimension(n_points), intent(out) :: included_n_reps
            !! Stores the count of non-NaN replicates (included ones)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
            !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)

        real(real64) :: bin_width, clamped_residual
        integer(int32) :: bin_idx, i_neighbor, i_rep, i_bin, included_reps, i_point
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
                    local(i_rep, clamped_residual, bin_idx) &
                    shared(filter_neighbors, n_reps, counts, neighborhood_residuals, shared_residual_range, bin_width) &
                    reduce(+:included_reps)
                ! Exclude neighbor if desired
                if (filter_neighbors) then
                    if (.not. neighbor_mask(i_neighbor, i_point)) cycle
                end if

                ! Count non-NaNs and assign the to a bin
                do i_rep = 1, n_reps
                    if (.not. ieee_is_nan(neighborhood_residuals(i_rep, i_neighbor, i_point))) then
                        ! clamp residual to histogram range
                        clamped_residual = clamp(neighborhood_residuals(i_rep, i_neighbor, i_point), min_val=-shared_residual_range, max_val=shared_residual_range)

                        ! assign bin to residual
                        bin_idx = min(n_bins, int( (clamped_residual + shared_residual_range) / bin_width ) + 1)
                        counts(i_point, bin_idx) = counts(i_point, bin_idx) + 1

                        included_reps = included_reps + 1
                    end if
                end do
            end do
            included_n_reps(i_point) = included_reps
        end do

        ! 2. calculate pmf
        do concurrent (i_bin = 1:n_bins)
            do concurrent (i_point = 1:n_points) shared(pmf, i_bin, included_n_reps, counts)
                if (included_n_reps(i_point) == 0) then
                    pmf(i_point, i_bin) = 0.0_real64
                else
                    pmf(i_point, i_bin) = real(counts(i_point, i_bin), real64) / real(included_n_reps(i_point), real64)
                end if
            end do
        end do
    end subroutine build_residual_histograms_helper

    !> Having the probabilities `pmf` from [[tox_data_integration(module):build_residual_histograms(interface)]], this subroutine computes the Jensen-Shannon divergence per reference point/neighbor
    pure module subroutine compute_divergence_per_reference_point(pmf_S1, pmf_S2, n_points, n_bins, js_divergences, ierr)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S1
            !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 1
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S2
            !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 2
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr)
        call validate_dimension_size(n_bins, ierr)
        call validate_all_in_range_real(pmf_S1, size(pmf_S1, kind=int32), ierr, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(pmf_S2, size(pmf_S2, kind=int32), ierr, min=0.0_real64, max=1.0_real64)

        if (is_err(ierr)) return

        call compute_divergence_per_reference_point_helper(pmf_S1, pmf_S2, n_points, n_bins, js_divergences)
    end subroutine compute_divergence_per_reference_point

    !> (no input validation) Having the probabilities `pmf` from [[tox_data_integration(module):build_residual_histograms(interface)]], this subroutine computes the Jensen-Shannon divergence per reference point/neighbor
    pure module subroutine compute_divergence_per_reference_point_helper(pmf_S1, pmf_S2, n_points, n_bins, js_divergences)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S1
            !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 1
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S2
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
                s1_val = pmf_S1(i_point, i_bin)
                s2_val = pmf_S2(i_point, i_bin)
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
            js_divergences(i_point) = 0.5_real64 * js_divergences(i_point)
        end do
    end subroutine compute_divergence_per_reference_point_helper

    !> Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
    pure module subroutine compute_weighted_global_divergence(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, ierr)
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

        call validate_dimension_size(n_points, ierr)
        call validate_all_in_range_real(js_divergences, size(js_divergences, kind=int32), ierr, min=0.0_real64)
        call validate_all_in_range_int(included_n_reps_S1, size(included_n_reps_S1, kind=int32), ierr, min=0_int32)
        call validate_all_in_range_int(included_n_reps_S2, size(included_n_reps_S2, kind=int32), ierr, min=0_int32)

        if (is_err(ierr)) return

        call compute_weighted_global_divergence_helper(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights)
    end subroutine compute_weighted_global_divergence

    !> (no input validation) Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
    pure module subroutine compute_weighted_global_divergence_helper(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights)
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

    !> Helper to run the pipeline `build_residual_histograms` \(\Rightarrow\) `compute_weighted_global_divergence`
    pure module subroutine jct_compute_jsd_pipeline_helper(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, neighbor_mask_S1, neighbor_mask_S2)
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
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S1
            !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S2
            !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

        call build_residual_histograms_helper(neighborhood_residuals_S1, n_reps_S1, n_neighbors, n_points, shared_residual_range, n_bins, tmp_counts, pmf_S1, included_n_reps_S1, neighbor_mask_S1)
        call build_residual_histograms_helper(neighborhood_residuals_S2, n_reps_S2, n_neighbors, n_points, shared_residual_range, n_bins, tmp_counts, pmf_S2, included_n_reps_S2, neighbor_mask_S2)
        call compute_divergence_per_reference_point_helper(pmf_S1, pmf_S2, n_points, n_bins, js_divergences)
        call compute_weighted_global_divergence_helper(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights)
    end subroutine jct_compute_jsd_pipeline_helper
end submodule tox_data_integration_jsd

!> C-compatible wrapper for [[tox_data_integration(module):determine_shared_residual_range(interface)]]
pure subroutine determine_shared_residual_range_expert_c( &
    abs_residual_pool, abs_residual_pool_perm, pool_size, &
    residual_range_quantile, shared_residual_range, ierr &
    ) bind(C, name="determine_shared_residual_range_expert_c")

    use tox_data_integration, only: determine_shared_residual_range
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: pool_size
        !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
    real(c_double), intent(in), target :: residual_range_quantile
        !! Quantile for determining the residual range, default: 95.0
    real(c_double), intent(out), target :: shared_residual_range
        !! Computed residual range (R)
    real(c_double), dimension(pool_size), intent(in), target :: abs_residual_pool
        !! The absolute residual values of the concatenated S1,S2 residuals
    integer(c_int), dimension(pool_size), intent(in), target :: abs_residual_pool_perm
        !! The permutation vector that sorts `abs_residual_pool`
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(pool_size)
    M_CHECK_NON_NULL(residual_range_quantile)
    M_CHECK_NON_NULL(shared_residual_range)
    M_CHECK_NON_NULL(abs_residual_pool)
    M_CHECK_NON_NULL(abs_residual_pool_perm)

    call determine_shared_residual_range( &
        abs_residual_pool, abs_residual_pool_perm, pool_size, &
        shared_residual_range, ierr, residual_range_quantile )

end subroutine determine_shared_residual_range_expert_c

!> C-compatible wrapper for [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
pure subroutine determine_shared_residual_range_c( &
    neighborhood_residuals_S1, neighborhood_residuals_S2, &
    n_reps_S1, n_reps_S2, n_neighbors, n_points, &
    residual_range_quantile, &
    shared_residual_range, ierr ) &
    bind(C, name="determine_shared_residual_range_c")

    use tox_data_integration, only: determine_shared_residual_range_alloc
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_reps_S1
        !! Number of replicates in study 1
    integer(c_int), intent(in), target :: n_reps_S2
        !! Number of replicates in study 2
    integer(c_int), intent(in), target :: n_neighbors
        !! Number of reference points (k)
    integer(c_int), intent(in), target :: n_points
        !! Number of reference points in the studies
    real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
        !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
        !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    real(c_double), intent(in), target :: residual_range_quantile
        !! Quantile for determining the residual range, default: 95.0
    real(c_double), intent(out), target :: shared_residual_range
        !! Computed residual range (R)
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_reps_S1)
    M_CHECK_NON_NULL(n_reps_S2)
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(neighborhood_residuals_S1)
    M_CHECK_NON_NULL(neighborhood_residuals_S2)
    M_CHECK_NON_NULL(residual_range_quantile)
    M_CHECK_NON_NULL(shared_residual_range)

    call determine_shared_residual_range_alloc( &
        neighborhood_residuals_S1, neighborhood_residuals_S2, &
        n_reps_S1, n_reps_S2, n_neighbors, n_points, &
        shared_residual_range, ierr, residual_range_quantile )

end subroutine determine_shared_residual_range_c

!> C-compatible wrapper for [[tox_data_integration(module):build_residual_histograms(interface)]]
pure subroutine build_residual_histograms_c( &
    neighborhood_residuals, &
    n_reps, n_neighbors, n_points, &
    shared_residual_range, &
    n_bins, &
    counts, pmf, included_n_reps, &
    ierr ) &
    bind(C, name="build_residual_histograms_c")

    use tox_data_integration, only: build_residual_histograms
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_reps
        !! Number of replicates in the study
    integer(c_int), intent(in), target :: n_neighbors
        !! Number of reference points (k)
    integer(c_int), intent(in), target :: n_points
        !! Number of reference points in the studies
    real(c_double), dimension(n_reps, n_neighbors, n_points), intent(in), target :: neighborhood_residuals
        !! Computed neighborhood residuals for a study ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    real(c_double), intent(in), target :: shared_residual_range
        !! Computed residual range (R) from [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
    integer(c_int), intent(in), target :: n_bins
        !! Number of equally sized histogram bins in range [-R,R]
    integer(c_int), dimension(n_points, n_bins), intent(out), target :: counts
        !! Absolute counts of a residual per bin
    real(c_double), dimension(n_points, n_bins), intent(out), target :: pmf
        !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
    integer(c_int), dimension(n_points), intent(out), target :: included_n_reps
        !! Stores the count of non-NaN replicates (included ones)
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(neighborhood_residuals)
    M_CHECK_NON_NULL(n_reps)
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(shared_residual_range)
    M_CHECK_NON_NULL(n_bins)
    M_CHECK_NON_NULL(counts)
    M_CHECK_NON_NULL(pmf)
    M_CHECK_NON_NULL(included_n_reps)

    call build_residual_histograms( &
        neighborhood_residuals, &
        n_reps, n_neighbors, n_points, &
        shared_residual_range, &
        n_bins, &
        counts, pmf, included_n_reps, &
        ierr )

end subroutine build_residual_histograms_c

!> C-compatible wrapper for [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
pure subroutine compute_divergence_per_reference_point_c( &
    pmf_S1, pmf_S2, &
    n_points, n_bins, &
    js_divergences, ierr ) &
    bind(C, name="compute_divergence_per_reference_point_c")

    use tox_data_integration, only: compute_divergence_per_reference_point
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_points
        !! Number of reference points (k)
    integer(c_int), intent(in), target :: n_bins
        !! Number of equally sized histogram bins in range [-R,R]
    real(c_double), dimension(n_points, n_bins), intent(in), target :: pmf_S1
        !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 1
    real(c_double), dimension(n_points, n_bins), intent(in), target :: pmf_S2
        !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 2
    real(c_double), dimension(n_points), intent(out), target :: js_divergences
        !! Jensen-Shannon divergence per reference point
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(n_bins)
    M_CHECK_NON_NULL(pmf_S1)
    M_CHECK_NON_NULL(pmf_S2)
    M_CHECK_NON_NULL(js_divergences)

    call compute_divergence_per_reference_point( &
        pmf_S1, pmf_S2, &
        n_points, n_bins, &
        js_divergences, ierr )

end subroutine compute_divergence_per_reference_point_c

!> C-compatible wrapper for [[tox_data_integration(module):compute_weighted_global_divergence(interface)]]
pure subroutine compute_weighted_global_divergence_c( &
    js_divergences, &
    n_points, &
    included_n_reps_S1, included_n_reps_S2, &
    global_js_divergence, weights, &
    ierr ) &
    bind(C, name="compute_weighted_global_divergence_c")

    use tox_data_integration, only: compute_weighted_global_divergence
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_points
        !! Number of reference points (k)
    real(c_double), dimension(n_points), intent(in), target :: js_divergences
        !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
    integer(c_int), dimension(n_points), intent(in), target :: included_n_reps_S1
        !! Count of non-NaN residuals (included ones) in study 1
    integer(c_int), dimension(n_points), intent(in), target :: included_n_reps_S2
        !! Count of non-NaN residuals (included ones) in study 2
    real(c_double), intent(out), target :: global_js_divergence
        !! Weighted global Jensen-Shannon divergence
    real(c_double), dimension(n_points), intent(out), target :: weights
        !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(js_divergences)
    M_CHECK_NON_NULL(included_n_reps_S1)
    M_CHECK_NON_NULL(included_n_reps_S2)
    M_CHECK_NON_NULL(global_js_divergence)
    M_CHECK_NON_NULL(weights)

    call compute_weighted_global_divergence( &
        js_divergences, n_points, &
        included_n_reps_S1, included_n_reps_S2, &
        global_js_divergence, weights, ierr )

end subroutine compute_weighted_global_divergence_c
