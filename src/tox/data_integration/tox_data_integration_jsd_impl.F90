#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation
!|
!| Implementations for the pipeline that turns neighborhood residuals -- as produced by the
!| preprocessing implementations -- into a JSD value. The generator turns the `*_impl` procedures
!| into the validating and allocating wrappers in module tox_data_integration_jsd;
!| `jct_compute_jsd_pipeline_helper` is an internal helper the permutation test and the
!| per-family analysis both drive, and carries no `_impl` suffix.
module tox_data_integration_jsd_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use f42_math_impl, only: clamp, is_close
    use f42_sort_impl, only: sort_array_heapsort
    use f42_utils, only: calc_percentile_helper
    M_IMPLICIT_NONE
contains

    !> summary: Compute the shared residual range [-R, R] from a pooled set of absolute residuals
    !| AUTHOR_FRANZ_ERIC_SILL
    !| This takes the pool already built; `determine_study_shared_residual_range` builds it from
    !| the neighborhood residuals of two studies first, if that is what is at hand.
    pure subroutine determine_shared_residual_range_impl(abs_residual_pool, abs_residual_pool_perm, pool_size, shared_residual_range, residual_range_quantile)
        integer(int32), intent(in) :: pool_size
            !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(95.0)
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), dimension(pool_size), intent(in) :: abs_residual_pool
            !! The absolute residual values of the concatenated S1,S2 residuals
            !! DM_ALLOW_NAN
        integer(int32), dimension(pool_size), intent(in) :: abs_residual_pool_perm
            !! The permutation vector that sorts `abs_residual_pool`
            !! DM_MIN(1_int32)
            !! DM_MAX(pool_size)

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

        call calc_percentile_helper(abs_residual_pool, abs_residual_pool_perm(:last_non_nan), actual_quantile, shared_residual_range)
    end subroutine determine_shared_residual_range_impl

    !> summary: Compute the shared residual range [-R, R] from the neighborhood residuals of two studies
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Pools the absolute residuals of both studies, sorts them, and takes the quantile exactly
    !| as `determine_shared_residual_range` does.
    pure subroutine determine_study_shared_residual_range_impl(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, &
                                                                 tmp_abs_residual_pool, tmp_abs_residual_pool_perm, shared_residual_range, residual_range_quantile)
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
        real(real64), dimension((n_reps_S1 + n_reps_S2)*n_neighbors*n_points), intent(out) :: tmp_abs_residual_pool
            !! Work array holding the pooled absolute residuals of both studies
        integer(int32), dimension((n_reps_S1 + n_reps_S2)*n_neighbors*n_points), intent(out) :: tmp_abs_residual_pool_perm
            !! Work array for the permutation that sorts `tmp_abs_residual_pool`
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(95.0)
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)

        integer(int32) :: i_rep, i_neighbor, pool_size, pool_idx, n_predecessors, i_point

        pool_size = (n_reps_S1 + n_reps_S2)*n_neighbors*n_points

        ! Collect the absolute residual values, laid out exactly as a
        ! (n_reps_S1 + n_reps_S2, n_neighbors, n_points) array would be
        do concurrent(i_point=1:n_points)
            do concurrent(i_neighbor=1:n_neighbors) local(n_predecessors) shared(i_point, n_reps_S1, n_reps_S2)
                n_predecessors = ((i_point - 1)*n_neighbors + (i_neighbor - 1))*(n_reps_S1 + n_reps_S2)
                do concurrent(i_rep=1:n_reps_S1) shared(n_predecessors, tmp_abs_residual_pool, tmp_abs_residual_pool_perm, neighborhood_residuals_S1, i_point, i_neighbor)
                    tmp_abs_residual_pool(n_predecessors + i_rep) = abs(neighborhood_residuals_S1(i_rep, i_neighbor, i_point))

                    tmp_abs_residual_pool_perm(n_predecessors + i_rep) = n_predecessors + i_rep
                end do

                do concurrent(i_rep=1:n_reps_S2) local(pool_idx) shared(n_predecessors, tmp_abs_residual_pool, tmp_abs_residual_pool_perm, neighborhood_residuals_S2, i_point, i_neighbor)
                    pool_idx = i_rep + n_reps_S1

                    tmp_abs_residual_pool(n_predecessors + pool_idx) = abs(neighborhood_residuals_S2(i_rep, i_neighbor, i_point))

                    tmp_abs_residual_pool_perm(n_predecessors + pool_idx) = n_predecessors + pool_idx
                end do
            end do
        end do

        call sort_array_heapsort(tmp_abs_residual_pool, tmp_abs_residual_pool_perm)

        call determine_shared_residual_range_impl(tmp_abs_residual_pool, tmp_abs_residual_pool_perm, pool_size, shared_residual_range, residual_range_quantile)
    end subroutine determine_study_shared_residual_range_impl

    !> summary: Summarize the neighborhood residuals in absolute histogram counts and probability mass functions
    !| AUTHOR_FRANZ_ERIC_SILL
    !| The probability mass function `pmf(residual, bin)` is actually a matrix.
    pure subroutine build_residual_histograms_impl(neighborhood_residuals, n_reps, n_neighbors, n_points, shared_residual_range, n_bins, counts, pmf, included_n_reps, neighbor_mask)
        integer(int32), intent(in) :: n_reps
            !! Number of replicates of the study
        integer(int32), intent(in) :: n_neighbors
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the study
        real(real64), dimension(n_reps, n_neighbors, n_points), intent(in) :: neighborhood_residuals
            !! Computed neighborhood residuals for a study, NaN is explicitly allowed for missing values
            !! DM_ALLOW_NAN
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
            !! DM_MIN(0.0_real64)
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

        ! Guard against a zero range (e.g. `determine_shared_residual_range_impl` returns 0.0 when all
        ! pooled residuals are NaN): fall back to a fixed bin width so every (clamped-to-zero) residual
        ! deterministically lands in a single bin instead of dividing by zero below.
        if (shared_residual_range <= 0.0_real64) then
            bin_width = 1.0_real64
        else
            bin_width = 2.0_real64*shared_residual_range/real(n_bins, real64)
        end if
        counts = 0_int32
        pmf = 0.0_real64

        filter_neighbors = present(neighbor_mask)

        ! 1. assign the bins to the residuals (increase the respective count)
        ! outer loop cannot be concurrent, as counts of same residuals and bins but different neighbors might be changed at the same time
        ! inner (i_neighbor) loop must be a plain sequential `do`: different neighbors can hit the same
        ! (i_point, bin_idx) count simultaneously, so `counts(i_point, bin_idx) = counts(i_point, bin_idx) + 1`
        ! would be a data race under `do concurrent`.
        do concurrent(i_point=1:n_points) local(included_reps) shared(included_n_reps)
            included_reps = 0_int32
            do i_neighbor = 1, n_neighbors
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
                        bin_idx = min(n_bins, int((clamped_residual + shared_residual_range)/bin_width) + 1)
                        counts(i_point, bin_idx) = counts(i_point, bin_idx) + 1

                        included_reps = included_reps + 1
                    end if
                end do
            end do
            included_n_reps(i_point) = included_reps
        end do

        ! 2. calculate pmf
        do concurrent(i_bin=1:n_bins)
            do concurrent(i_point=1:n_points) shared(pmf, i_bin, included_n_reps, counts)
                if (included_n_reps(i_point) == 0) then
                    pmf(i_point, i_bin) = 0.0_real64
                else
                    pmf(i_point, i_bin) = real(counts(i_point, i_bin), real64)/real(included_n_reps(i_point), real64)
                end if
            end do
        end do
    end subroutine build_residual_histograms_impl

    !> summary: Compute the Jensen-Shannon divergence per reference point from two histograms
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Takes the probabilities `pmf` produced by `build_residual_histograms`.
    pure subroutine compute_divergence_per_reference_point_impl(pmf_S1, pmf_S2, n_points, n_bins, js_divergences)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S1
            !! Computed normalized histogram counts for study 1
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S2
            !! Computed normalized histogram counts for study 2
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
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
            do concurrent(i_point=1:n_points) local(s1_val, s2_val, S_mean) shared(i_bin, pmf_S1, pmf_S2, js_divergences)
                s1_val = pmf_S1(i_point, i_bin)
                s2_val = pmf_S2(i_point, i_bin)
                S_mean = 0.5_real64*(s1_val + s2_val)

                if (.not. is_close(S_mean, 0.0_real64)) then
                    if (s1_val > 0.0_real64) then
                        js_divergences(i_point) = js_divergences(i_point) + s1_val*log(s1_val/S_mean)
                    end if

                    if (s2_val > 0.0_real64) then
                        js_divergences(i_point) = js_divergences(i_point) + s2_val*log(s2_val/S_mean)
                    end if
                end if
            end do
        end do

        ! 2. Compute the js_divergences
        do concurrent(i_point=1:n_points) shared(js_divergences)
            js_divergences(i_point) = 0.5_real64*js_divergences(i_point)
        end do
    end subroutine compute_divergence_per_reference_point_impl

    !> summary: Compute the global weighted Jensen-Shannon divergence from the per-neighbor divergences
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Takes the divergences produced by `compute_divergence_per_reference_point`.
    pure subroutine compute_weighted_global_divergence_impl(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        real(real64), dimension(n_points), intent(in) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
            !! DM_MIN(0.0_real64)
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
            !! DM_MIN(0_int32)
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`

        integer(int32) :: i_point, included_reps
        real(real64) :: total_sample_count

        global_js_divergence = 0.0_real64

        total_sample_count = real(sum(included_n_reps_S1) + sum(included_n_reps_S2), real64)

        ! Guard against a zero total sample count (no replicate had a non-NaN residual at any reference
        ! point in either study): fall back to all-zero weights instead of dividing by zero below.
        if (is_close(total_sample_count, 0.0_real64)) then
            weights = 0.0_real64
        else
            ! Calculate the global Jensen-Shannon divergence
            do concurrent(i_point=1:n_points) local(included_reps) shared(weights, included_n_reps_S1, included_n_reps_S2, total_sample_count) reduce(+:global_js_divergence)
                included_reps = included_n_reps_S1(i_point) + included_n_reps_S2(i_point)

                weights(i_point) = real(included_reps, real64)/total_sample_count

                global_js_divergence = global_js_divergence + weights(i_point)*js_divergences(i_point)
            end do
        end if
    end subroutine compute_weighted_global_divergence_impl

    !> summary: Run the pipeline build_residual_histograms => compute_weighted_global_divergence
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Internal helper: the permutation test and the per-family analysis both drive this.
    pure subroutine jct_compute_jsd_pipeline_helper(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, neighbor_mask_S1, neighbor_mask_S2)
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
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S1
            !! Normalized histogram counts for study 1
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S2
            !! Normalized histogram counts for study 2
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array for the histogram counts
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

        call build_residual_histograms_impl(neighborhood_residuals_S1, n_reps_S1, n_neighbors, n_points, shared_residual_range, n_bins, tmp_counts, pmf_S1, included_n_reps_S1, neighbor_mask_S1)
        call build_residual_histograms_impl(neighborhood_residuals_S2, n_reps_S2, n_neighbors, n_points, shared_residual_range, n_bins, tmp_counts, pmf_S2, included_n_reps_S2, neighbor_mask_S2)
        call compute_divergence_per_reference_point_impl(pmf_S1, pmf_S2, n_points, n_bins, js_divergences)
        call compute_weighted_global_divergence_impl(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights)
    end subroutine jct_compute_jsd_pipeline_helper
end module tox_data_integration_jsd_impl
